#!/usr/bin/env python
# coding: utf-8

#### Dependencies
import argparse
import os
import re
import requests
import pandas as pd
from Bio import PDB

#### Config
BASE_DIR = "data/pdbs_query"
CIF_DIR = os.path.join(BASE_DIR, "structures")
os.makedirs(BASE_DIR, exist_ok=True)
os.makedirs(CIF_DIR, exist_ok=True)

#### Functions

def uniprot_to_pdb_ids(uniprot_id):
    """Query RCSB Search API for PDB IDs associated with a UniProt accession.

    Iterates through paginated results until all matching PDB IDs are retrieved.

    Parameters
    ----------
    uniprot_id : str
        UniProt accession (e.g. "Q9UNS1").

    Returns
    -------
    list of str
        PDB entry identifiers linked to the UniProt ID.
    """
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    all_pdb_ids = []
    rows_per_page = 10000
    start = 0

    while True:
        query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                    "operator": "exact_match",
                    "value": uniprot_id
                }
            },
            "request_options": {
                "paginate": {
                    "start": start,
                    "rows": rows_per_page
                }
            },
            "return_type": "entry"
        }

        try:
            response = requests.post(
                url,
                json=query,
                headers={"Content-Type": "application/json"},
                timeout=60
            )
            response.raise_for_status()
            if not response.text:
                break
            data = response.json()

            total_count = data.get("total_count", 0)
            result_set = data.get("result_set", [])
            page_ids = [entry["identifier"] for entry in result_set]
            all_pdb_ids.extend(page_ids)

            if len(all_pdb_ids) >= total_count or len(page_ids) == 0:
                break

            start += rows_per_page

        except Exception as e:
            print(f"Warning: Request failed at offset {start}: {e}")
            break

    return all_pdb_ids


def get_pdb_structures_info(pdb_ids):
    """Fetch metadata for a list of PDB entries via the RCSB Data API.

    Parameters
    ----------
    pdb_ids : iterable of str
        PDB entry identifiers.

    Returns
    -------
    pandas.DataFrame
        Columns: pdb_id, resolution, residue_coverage, num_chains, method.
    """
    rows = []
    for pdb_id in pdb_ids:
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            data = response.json()

            resolution = data.get("rcsb_entry_info", {}).get("resolution_combined")
            if isinstance(resolution, list):
                resolution = resolution[0] if resolution else None

            method = None
            exptl = data.get("exptl", [])
            if exptl:
                method = exptl[0].get("method")

            num_chains = data.get(
                "rcsb_entry_info", {}
            ).get("polymer_entity_instance_count")
            residue_coverage = data.get(
                "rcsb_entry_info", {}
            ).get("deposited_polymer_monomer_count")

            rows.append({
                "pdb_id": pdb_id,
                "resolution": resolution,
                "residue_coverage": residue_coverage,
                "num_chains": num_chains,
                "method": method
            })
        except Exception:
            pass

    return pd.DataFrame(
        rows,
        columns=["pdb_id", "resolution", "residue_coverage", "num_chains", "method"]
    )


def download_cif(pdb_id, out_dir):
    """Download a mmCIF file from RCSB and save it to disk.

    Skips download if the file already exists.

    Parameters
    ----------
    pdb_id : str
        PDB entry identifier.
    out_dir : str
        Directory where the .cif file will be written.

    Returns
    -------
    str or None
        Path to the saved file, or None if the download failed.
    """
    path = os.path.join(out_dir, f"{pdb_id}.cif")
    if os.path.exists(path):
        return path

    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    try:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        with open(path, "w") as fh:
            fh.write(response.text)
        return path
    except Exception:
        return None


def parse_dbref_records(file_path):
    """Extract DBREF records from a PDB-format file.

    Parameters
    ----------
    file_path : str
        Path to a PDB file.

    Returns
    -------
    list of dict
        Each dict contains ``chain_id`` and ``uniprot_acc`` keys.
    """
    dbref_data = []
    try:
        with open(file_path, "r") as fh:
            for line in fh:
                if line.startswith("DBREF"):
                    chain_id = line[12].strip() if len(line) > 12 else ""
                    database = line[26:32].strip() if len(line) > 32 else ""
                    db_accession = line[33:41].strip() if len(line) > 41 else ""
                    if database == "UNP" and db_accession:
                        dbref_data.append({
                            "chain_id": chain_id,
                            "uniprot_acc": db_accession
                        })
    except Exception:
        pass

    return dbref_data


def get_chain_identity(cif_file_path, uniprot_id):
    """Parse a mmCIF file and classify each chain relative to a target UniProt ID.

    Chains are labelled target when their DBREF UniProt matches
    uniprot_id, other when a different UniProt is mapped, and
    unknown when no mapping is found.

    UniProt-numbered ranges come from _struct_ref_seq.db_align_beg /
    db_align_end (NOT from seq_align_beg/end, which are entity/PDB-sequence
    numbering, and NOT from a "_struct_ref_seq_seq" category, which does not
    exist in the mmCIF dictionary).

    Returns
    -------
    pandas.DataFrame
        Columns include uniprot_id, pdb_id, chain_id, chain_uniprot,
        identity, uniprot_ranges, pdb_auth_ranges, num_residues, num_atoms,
        resolution, method, title.
    """
    parser = PDB.MMCIFParser(QUIET=True)
    pdb_id = os.path.basename(cif_file_path).replace(".cif", "")

    mmcif_dict = PDB.MMCIF2Dict.MMCIF2Dict(cif_file_path)

    def as_list(x):
        if x is None:
            return []
        return x if isinstance(x, list) else [x]

    strand_ids   = as_list(mmcif_dict.get("_struct_ref_seq.pdbx_strand_id"))
    db_accessions = as_list(mmcif_dict.get("_struct_ref_seq.pdbx_db_accession"))
    db_align_begs = as_list(mmcif_dict.get("_struct_ref_seq.db_align_beg"))
    db_align_ends = as_list(mmcif_dict.get("_struct_ref_seq.db_align_end"))
    auth_begs     = as_list(mmcif_dict.get("_struct_ref_seq.pdbx_auth_seq_align_beg"))
    auth_ends     = as_list(mmcif_dict.get("_struct_ref_seq.pdbx_auth_seq_align_end"))

    n = len(strand_ids)
    # Pad any missing/mismatched columns so zip doesn't silently truncate.
    def pad(lst):
        return lst + [None] * (n - len(lst)) if len(lst) < n else lst

    db_accessions = pad(db_accessions)
    db_align_begs = pad(db_align_begs)
    db_align_ends = pad(db_align_ends)
    auth_begs = pad(auth_begs)
    auth_ends = pad(auth_ends)

    chain_to_uniprot = {}
    chain_to_uniprot_ranges = {}
    chain_to_auth_ranges = {}

    for strand_id, db_acc, ub, ue, ab, ae in zip(
        strand_ids, db_accessions, db_align_begs, db_align_ends, auth_begs, auth_ends
    ):
        if strand_id is None:
            continue
        for sid in strand_id.split(","):
            sid = sid.strip()
            if not sid:
                continue
            if db_acc:
                chain_to_uniprot[sid] = db_acc.strip()

            chain_to_uniprot_ranges.setdefault(sid, [])
            chain_to_auth_ranges.setdefault(sid, [])
            try:
                chain_to_uniprot_ranges[sid].append((int(ub), int(ue)))
            except (TypeError, ValueError):
                pass
            try:
                chain_to_auth_ranges[sid].append((int(ab), int(ae)))
            except (TypeError, ValueError):
                pass

    def format_ranges(ranges):
        if not ranges:
            return None
        sorted_ranges = sorted(ranges, key=lambda x: x[0])
        merged = [list(sorted_ranges[0])]
        for b, e in sorted_ranges[1:]:
            if b <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([b, e])
        return ";".join(f"{m[0]}-{m[1]}" for m in merged)

    structure = parser.get_structure(pdb_id, cif_file_path)
    header = structure.header
    resolution = header.get("resolution", None)
    method = header.get("structure_method", None)

    # Retrieve the descriptive title from the mmCIF dictionary
    title = mmcif_dict.get("_struct.title", None)
    if isinstance(title, list):
        title = title[0] if title else None
    elif isinstance(title, str):
        title = title.strip() or None
    else:
        title = None

    rows = []
    for model in structure:
        for chain in model:
            chain_id = chain.id
            residues = list(chain.get_residues())
            num_residues = len(residues)
            num_atoms = sum(len(list(r.get_atoms())) for r in residues)

            chain_uniprot = chain_to_uniprot.get(chain_id, None)
            if chain_uniprot == uniprot_id:
                identity = "target"
            elif chain_uniprot:
                identity = "other"
            else:
                identity = "unknown"

            uniprot_ranges = format_ranges(chain_to_uniprot_ranges.get(chain_id, []))
            pdb_auth_ranges = format_ranges(chain_to_auth_ranges.get(chain_id, []))

            rows.append({
                "uniprot_id": uniprot_id,
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "chain_uniprot": chain_uniprot,
                "identity": identity,
                "uniprot_ranges": uniprot_ranges,
                "pdb_auth_ranges": pdb_auth_ranges,
                "num_residues": num_residues,
                "num_atoms": num_atoms,
                "resolution": resolution,
                "method": method,
                "title": title
            })
        break  # only first model needed for chain-level mapping

    return pd.DataFrame(rows)
    

#### Main
def main():
    parser = argparse.ArgumentParser(
        description="Map UniProt ID to PDB structures and chain identities"
    )
    parser.add_argument(
        "--uniprot_id",
        required=True,
        help="UniProt accession ID"
    )
    args = parser.parse_args()

    uniprot_id = args.uniprot_id
    pdbs_file = os.path.join(BASE_DIR, f"{uniprot_id}_pdbs.tsv")
    chains_file = os.path.join(BASE_DIR, f"{uniprot_id}_chains.tsv")

    #### 1. Get PDB IDs
    if os.path.exists(pdbs_file):
        print(f"Loading existing PDB list from {pdbs_file}")
        df_pdbs = pd.read_csv(pdbs_file, sep="\t")
    else:
        print(f"Querying PDB IDs for {uniprot_id}...")
        pdb_ids = uniprot_to_pdb_ids(uniprot_id)
        if not pdb_ids:
            print(f"No PDB IDs found for {uniprot_id}")
            exit(0)
        print(f"Found {len(pdb_ids)} PDB IDs")
        df_pdbs = get_pdb_structures_info(pdb_ids)
        df_pdbs["uniprot_id"] = uniprot_id
        df_pdbs = df_pdbs[[
            "uniprot_id", "pdb_id", "resolution", "residue_coverage",
            "num_chains", "method"
        ]]
        df_pdbs.to_csv(pdbs_file, sep="\t", index=False)
        print(f"Saved PDB metadata to {pdbs_file}")

    #### 2. Download CIFs and extract chain identity data
    if os.path.exists(chains_file):
        print(f"Loading existing chain data from {chains_file}")
        df_chains = pd.read_csv(chains_file, sep="\t")
    else:
        print("Downloading CIF structures and parsing chain data...")
        all_chain_rows = []
        for pdb_id in df_pdbs["pdb_id"].unique():
            cif_path = download_cif(pdb_id, CIF_DIR)
            if cif_path:
                chain_df = get_chain_identity(cif_path, uniprot_id)
                all_chain_rows.append(chain_df)
        if all_chain_rows:
            df_chains = pd.concat(all_chain_rows, ignore_index=True)
        else:
            df_chains = pd.DataFrame(columns=[
                "uniprot_id", "pdb_id", "chain_id", "chain_uniprot",
                "identity", "pdb_to_uniprot_ranges", "num_residues",
                "num_atoms", "resolution", "method"
            ])
        df_chains.to_csv(chains_file, sep="\t", index=False)
        print(f"Saved chain identity data to {chains_file}")

    #### Output
    print("\n--- PDB Metadata ---")
    print(df_pdbs.head())
    print(f"\nTotal PDB entries: {len(df_pdbs)}")

    print("\n--- Chain Identity Data ---")
    print(df_chains.head())
    print(f"\nTotal chain records: {len(df_chains)}")
    print(f"\nColumns: {list(df_chains.columns)}")


if __name__ == "__main__":
    main()

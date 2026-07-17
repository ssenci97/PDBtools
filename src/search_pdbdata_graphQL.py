#!/usr/bin/env python3
# coding: utf-8

import argparse
import os
import requests
import pandas as pd


# ─── TUNABLES ────────────────────────────────────────────────────────────────
ROWS_PER_PAGE = 10000
GRAPHQL_BATCH_SIZE = 100
GRAPHQL_TIMEOUT = 60
SEARCH_TIMEOUT = 60
# ─────────────────────────────────────────────────────────────────────────────


def search_rcsb(query_text, min_year):
    """Query RCSB Search API v2 combining a free-text search with a release year limit."""
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    all_pdb_ids = []
    start = 0

    search_payload = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {"value": query_text}
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_accession_info.initial_release_date",
                        "operator": "greater_or_equal",
                        "value": f"{min_year}-01-01"
                    }
                }
            ]
        },
        "request_options": {
            "paginate": {
                "start": start,
                "rows": ROWS_PER_PAGE
            }
        },
        "return_type": "entry"
    }

    while True:
        try:
            response = requests.post(
                url,
                json=search_payload,
                headers={"Content-Type": "application/json"},
                timeout=SEARCH_TIMEOUT
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

            start += ROWS_PER_PAGE
            search_payload["request_options"]["paginate"]["start"] = start

        except Exception as e:
            print(f"Warning: Search failed at offset {start}: {e}")
            break

    return all_pdb_ids


def fetch_metadata(pdb_ids):
    """Fetch entry-level metadata + stoichiometry + symmetry via GraphQL."""
    if not pdb_ids:
        return pd.DataFrame()

    url = "https://data.rcsb.org/graphql"

    graphql_query = """
    query GetPdbDetails($ids: [String!]!) {
      entries(entry_ids: $ids) {
        rcsb_id
        rcsb_accession_info {
          initial_release_date
        }
        rcsb_entry_info {
          resolution_combined
          polymer_entity_count
          deposited_polymer_entity_instance_count
          polymer_entity_count_protein
          polymer_entity_count_DNA
          polymer_entity_count_RNA
        }
        exptl {
          method
        }
        struct {
          title
        }
        polymer_entities {
          rcsb_id
          rcsb_polymer_entity_container_identifiers {
            auth_asym_ids
          }
          entity_poly {
            pdbx_strand_id
            rcsb_entity_polymer_type
          }
        }
        assemblies {
          rcsb_id
          rcsb_assembly_info {
            assembly_id
            polymer_entity_count
            polymer_entity_instance_count
            polymer_entity_count_protein
            polymer_entity_count_DNA
            polymer_entity_count_RNA
            polymer_composition
            selected_polymer_entity_types
          }
          rcsb_struct_symmetry {
            symbol
            type
            oligomeric_state
            stoichiometry
          }
        }
      }
    }
    """

    all_rows = []

    for i in range(0, len(pdb_ids), GRAPHQL_BATCH_SIZE):
        batch = pdb_ids[i:i + GRAPHQL_BATCH_SIZE]
        payload = {
            "query": graphql_query,
            "variables": {"ids": batch}
        }
        try:
            response = requests.post(url, json=payload, timeout=GRAPHQL_TIMEOUT)
            response.raise_for_status()
            entries = response.json().get("data", {}).get("entries") or []

            for entry in entries:
                pdb_id = entry.get("rcsb_id")
                release_date = entry.get("rcsb_accession_info", {}).get("initial_release_date")
                resolutions = entry.get("rcsb_entry_info", {}).get("resolution_combined")
                resolution = resolutions[0] if resolutions else None
                exptl = entry.get("exptl") or []
                method = exptl[0].get("method") if exptl else None
                title = entry.get("struct", {}).get("title")

                info = entry.get("rcsb_entry_info") or {}
                polymer_entities = entry.get("polymer_entities") or []

                # Per-entity chain info
                entity_summaries = []
                for ent in polymer_entities:
                    cid = ent.get("rcsb_id")
                    chain_ids = ent.get("rcsb_polymer_entity_container_identifiers", {}).get("auth_asym_ids") or []
                    strand_id = ent.get("entity_poly", {}).get("pdbx_strand_id")
                    ptype = ent.get("entity_poly", {}).get("rcsb_entity_polymer_type")
                    entity_summaries.append(
                        f"{cid}({ptype or '?'})={','.join(chain_ids) or strand_id or '?'}"
                    )

                # Assembly / stoichiometry / symmetry info
                assemblies = entry.get("assemblies") or []
                assembly_summaries = []
                symmetry_summaries = []
                stoichiometry_summaries = []

                for asm in assemblies:
                    asm_id = asm.get("rcsb_id")
                    ainfo = asm.get("rcsb_assembly_info") or {}
                    sym_list = asm.get("rcsb_struct_symmetry") or []

                    # Assembly composition summary
                    asm_summary = (
                        f"{asm_id}: "
                        f"entities={ainfo.get('polymer_entity_count')}, "
                        f"instances={ainfo.get('polymer_entity_instance_count')}, "
                        f"composition={ainfo.get('polymer_composition') or '?'}"
                    )
                    assembly_summaries.append(asm_summary)

                    # Symmetry info
                    for sym in sym_list:
                        sym_summary = (
                            f"{sym.get('symbol')} ({sym.get('type')}) — "
                            f"{sym.get('oligomeric_state')}"
                        )
                        symmetry_summaries.append(sym_summary)

                        # Stoichiometry (e.g., ["A2", "B2"] or ["A4"])
                        stoich = sym.get("stoichiometry") or []
                        if stoich:
                            stoichiometry_summaries.append(",".join(stoich))

                all_rows.append({
                    "pdb_id": pdb_id,
                    "release_date": release_date,
                    "resolution": resolution,
                    "method": method,
                    "title": title,
                    # Deposited model counts
                    "polymer_entity_count": info.get("polymer_entity_count"),
                    "polymer_instance_count": info.get("deposited_polymer_entity_instance_count"),
                    "protein_entity_count": info.get("polymer_entity_count_protein"),
                    "dna_entity_count": info.get("polymer_entity_count_DNA"),
                    "rna_entity_count": info.get("polymer_entity_count_RNA"),
                    "entity_chain_map": "; ".join(entity_summaries) if entity_summaries else None,
                    # Assembly-level info
                    "assembly_details": "; ".join(assembly_summaries) if assembly_summaries else None,
                    "global_symmetry": "; ".join(symmetry_summaries) if symmetry_summaries else None,
                    "global_stoichiometry": "; ".join(stoichiometry_summaries) if stoichiometry_summaries else None,
                })
        except Exception as e:
            print(f"Warning: Metadata fetch failed for batch {i}: {e}")

    return pd.DataFrame(all_rows)


def main():
    parser = argparse.ArgumentParser(
        description="Search RCSB with full-text queries and fetch entry metadata + stoichiometry + symmetry."
    )
    parser.add_argument("--query", required=True, help="Search query (e.g. 'kinase', 'ribosome')")
    parser.add_argument("--min_year", type=int, default=2023, help="Keep PDBs released from this year onward")
    parser.add_argument("--outdir", default="data/pdbs_query", help="Output directory")
    args = parser.parse_args()

    safe_query_name = "".join(c if c.isalnum() else "_" for c in args.query).strip("_")
    os.makedirs(args.outdir, exist_ok=True)
    out_file = os.path.join(args.outdir, f"{safe_query_name}_pdbs.tsv")

    print(f"Searching RCSB for: '{args.query}' (Released {args.min_year}+)...")
    pdb_ids = search_rcsb(args.query, args.min_year)
    if not pdb_ids:
        print("No matching structures found.")
        return
    print(f"Found {len(pdb_ids)} structures.")

    print("Fetching metadata + stoichiometry + symmetry via GraphQL...")
    df = fetch_metadata(pdb_ids)

    if not df.empty:
        df.to_csv(out_file, sep="\t", index=False)
        print(f"Saved PDB metadata to {out_file}")

        print(f"\n--- Preview of Results ({len(df)} total) ---")
        preview_cols = [
            "pdb_id", "release_date", "resolution", "method",
            "polymer_entity_count", "polymer_instance_count",
            "global_stoichiometry", "global_symmetry"
        ]
        print(df[[c for c in preview_cols if c in df.columns]].head(10))
    else:
        print("Failed to compile metadata.")


if __name__ == "__main__":
    main()


#!/usr/bin/env python3
import sys
import argparse
from pathlib import Path
import urllib.request
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

def download_mmcif(pdb_id, outdir):
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    filepath = Path(outdir) / f"{pdb_id}.cif"
    urllib.request.urlretrieve(url, filepath)
    return filepath

def three_to_one(resname):
    try:
        return seq1(resname)
    except KeyError:
        return 'X'

def get_entity_id_from_chain(chain_id, mmcif_dict, fallback_entity=None):
    if '_struct_asym.id' in mmcif_dict:
        asym_ids = mmcif_dict['_struct_asym.id']
        entity_ids = mmcif_dict.get('_struct_asym.entity_id', [])
        auth_asym_ids = mmcif_dict.get('_struct_asym.pdbx_auth_asym_id', asym_ids)
        for i, asym_id in enumerate(asym_ids):
            if asym_id == chain_id or (i < len(auth_asym_ids) and auth_asym_ids[i] == chain_id):
                if i < len(entity_ids) and entity_ids[i] not in ('?', None, ''):
                    return entity_ids[i]
    if '_entity_poly.entity_id' in mmcif_dict:
        entity_ids = mmcif_dict['_entity_poly.entity_id']
        pdbx_strand_ids = mmcif_dict.get('_entity_poly.pdbx_strand_id', [])
        for i, entity_id in enumerate(entity_ids):
            if i < len(pdbx_strand_ids):
                chains_for_entity = pdbx_strand_ids[i].split(',')
                if chain_id in chains_for_entity:
                    return entity_id
    if '_struct_ref_seq.pdbx_strand_id' in mmcif_dict:
        strand_ids = mmcif_dict['_struct_ref_seq.pdbx_strand_id']
        entity_ids = mmcif_dict.get('_struct_ref_seq.entity_id', [])
        for i, strand_id in enumerate(strand_ids):
            if strand_id == chain_id and i < len(entity_ids) and entity_ids[i] not in ('?', None):
                return entity_ids[i]
    return fallback_entity or "?"

def parse_db_references(mmcif_dict):
    db_references = {}
    if '_struct_ref.id' in mmcif_dict:
        struct_ref_ids = mmcif_dict['_struct_ref.id']
        db_names = mmcif_dict.get('_struct_ref.db_name', [])
        db_codes = mmcif_dict.get('_struct_ref.pdbx_db_accession', [])
        db_accessions = mmcif_dict.get('_struct_ref.db_code', [])
        entity_ids = mmcif_dict.get('_struct_ref.entity_id', [])
        for i, ref_id in enumerate(struct_ref_ids):
            entity_id = entity_ids[i] if i < len(entity_ids) else None
            if not entity_id or entity_id == '?':
                if '_struct_ref_seq.ref_id' in mmcif_dict:
                    seq_ref_ids = mmcif_dict['_struct_ref_seq.ref_id']
                    seq_entities = mmcif_dict.get('_struct_ref_seq.entity_id', [])
                    for j, srid in enumerate(seq_ref_ids):
                        if srid == ref_id and j < len(seq_entities) and seq_entities[j] not in ('?', None):
                            entity_id = seq_entities[j]
                            break
            if entity_id and entity_id != '?':
                db_name = db_names[i] if i < len(db_names) else None
                db_code = db_codes[i] if i < len(db_codes) else None
                if not db_code and i < len(db_accessions):
                    db_code = db_accessions[i]
                if db_name and db_code:
                    if db_name.upper() in ['UNP', 'UNIPROT']:
                        db_name = 'UniProt'
                    if entity_id not in db_references:
                        db_references[entity_id] = []
                    ref_entry = {'db_name': db_name, 'db_code': db_code}
                    if ref_entry not in db_references[entity_id]:
                        db_references[entity_id].append(ref_entry)
    if '_database_2.database_id' in mmcif_dict:
        db_ids = mmcif_dict['_database_2.database_id']
        db_codes = mmcif_dict.get('_database_2.database_code', [])
        protein_entities = []
        if '_entity_poly.entity_id' in mmcif_dict:
            for eid in mmcif_dict['_entity_poly.entity_id']:
                protein_entities.append(eid)
        elif '_entity.id' in mmcif_dict:
            for i, eid in enumerate(mmcif_dict['_entity.id']):
                etype = mmcif_dict.get('_entity.type', [None]*len(mmcif_dict['_entity.id']))[i]
                if etype == 'polymer':
                    protein_entities.append(eid)
        for db_idx, entity_id in enumerate(protein_entities):
            if db_idx < len(db_ids):
                db_id = db_ids[db_idx]
                db_code = db_codes[db_idx] if db_idx < len(db_codes) else None
                if db_id.upper() in ['UNP', 'UNIPROT', 'SWS'] and db_code:
                    if entity_id not in db_references:
                        db_references[entity_id] = []
                    if not any(r['db_code'] == db_code for r in db_references[entity_id]):
                        db_references[entity_id].append({'db_name': 'UniProt', 'db_code': db_code})
    return db_references

def parse_with_biopython(filepath, pdb_id):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(pdb_id.upper(), str(filepath))
    mmcif_dict = MMCIF2Dict(str(filepath))
    entity_info = {}
    if '_entity.id' in mmcif_dict:
        for i, entity_id in enumerate(mmcif_dict['_entity.id']):
            entity_type = mmcif_dict.get('_entity.type', [None] * len(mmcif_dict['_entity.id']))[i]
            entity_info[entity_id] = {'type': entity_type}
    db_references = parse_db_references(mmcif_dict)
    chain_to_entity = {}
    auth_to_label_chain = {}
    if '_struct_asym.id' in mmcif_dict:
        for i, label_asym_id in enumerate(mmcif_dict['_struct_asym.id']):
            entity_id = mmcif_dict.get('_struct_asym.entity_id', [None] * len(mmcif_dict['_struct_asym.id']))[i]
            auth_asym_id = mmcif_dict.get('_struct_asym.pdbx_auth_asym_id', [label_asym_id])[i] if '_struct_asym.pdbx_auth_asym_id' in mmcif_dict else label_asym_id
            if entity_id and entity_id != '?':
                chain_to_entity[label_asym_id] = entity_id
                chain_to_entity[auth_asym_id] = entity_id
            auth_to_label_chain[auth_asym_id] = label_asym_id
    if not chain_to_entity and '_entity_poly.entity_id' in mmcif_dict:
        poly_entities = mmcif_dict['_entity_poly.entity_id']
        strand_ids = mmcif_dict.get('_entity_poly.pdbx_strand_id', [])
        for i, entity_id in enumerate(poly_entities):
            if i < len(strand_ids):
                for chain in strand_ids[i].split(','):
                    chain_to_entity[chain.strip()] = entity_id
    chains_info = {}
    for model in structure:
        for chain in model:
            chain_id = chain.id.strip()
            if not chain_id:
                continue
            entity_id = chain_to_entity.get(chain_id)
            if not entity_id or entity_id == '?':
                entity_id = get_entity_id_from_chain(chain_id, mmcif_dict)
            label_asym_id = auth_to_label_chain.get(chain_id, chain_id)
            db_refs = db_references.get(entity_id, [])
            if not db_refs and '_struct_ref_seq.pdbx_strand_id' in mmcif_dict:
                seq_strands = mmcif_dict['_struct_ref_seq.pdbx_strand_id']
                seq_entities = mmcif_dict.get('_struct_ref_seq.entity_id', [])
                for i, strand in enumerate(seq_strands):
                    if strand == chain_id and i < len(seq_entities) and seq_entities[i] not in ('?', None):
                        entity_id = seq_entities[i]
                        db_refs = db_references.get(entity_id, [])
                        break
            entity_type = entity_info.get(entity_id, {}).get('type', '?')
            residues = []
            ca_coords = []
            auth_seqids = []
            all_atoms = []
            for residue in chain:
                resname = residue.resname.strip()
                if not resname or resname == 'HOH':
                    continue
                hetflag = residue.id[0].strip()
                seq_num = residue.id[1]
                icode = residue.id[2].strip()
                auth_seq = str(seq_num) + icode
                if entity_type == 'polymer':
                    aa = three_to_one(resname)
                else:
                    aa = resname
                for atom in residue:
                    atom_name = atom.name.strip()
                    x, y, z = atom.coord
                    occupancy = atom.get_occupancy()
                    b_factor = atom.get_bfactor()
                    element = atom.element.strip() if atom.element else ''
                    all_atoms.append({
                        'atom_name': atom_name,
                        'resname': resname,
                        'hetflag': hetflag,
                        'seq_num': seq_num,
                        'icode': icode,
                        'auth_seq': auth_seq,
                        'x': x, 'y': y, 'z': z,
                        'occupancy': occupancy,
                        'b_factor': b_factor,
                        'element': element,
                    })
                ca_x = ca_y = ca_z = None
                if 'CA' in residue:
                    atom = residue['CA']
                    ca_x, ca_y, ca_z = atom.coord
                residues.append(aa)
                if ca_x is not None:
                    ca_coords.append((auth_seq, ca_x, ca_y, ca_z))
                auth_seqids.append(auth_seq)
            if not residues:
                continue
            gapped = []
            has_gaps = False
            prev_num = None
            if entity_type == 'polymer':
                for i, seqid in enumerate(auth_seqids):
                    aa = residues[i]
                    try:
                        num = int(''.join(c for c in seqid if c.isdigit()))
                    except ValueError:
                        num = None
                    if prev_num is not None and num is not None:
                        gap = num - prev_num - 1
                        if gap > 0:
                            gapped.extend(['X'] * gap)
                            has_gaps = True
                    gapped.append(aa)
                    if num is not None:
                        prev_num = num
                gapped_seq = ''.join(gapped)
            else:
                gapped_seq = ''
                has_gaps = False
            plain_seq = ''.join(residues)
            chains_info[chain_id] = {
                'auth_asym_id': chain_id,
                'label_asym_id': label_asym_id,
                'entity_id': entity_id,
                'entity_type': entity_type,
                'db_references': db_refs,
                'sequence': plain_seq,
                'gapped_sequence': gapped_seq,
                'has_gaps': has_gaps,
                'ca_positions': ca_coords,
                'auth_seqids': auth_seqids,
                'all_atoms': all_atoms,
            }
    return chains_info

def write_fasta_files(chains_info, pdb_id, outdir):
    seen_sequences = set()
    written_count = 0
    for chain_id, info in sorted(chains_info.items()):
        if info['has_gaps'] and info['gapped_sequence']:
            sequence = info['gapped_sequence']
            gaps_flag = 'yes'
        elif info['sequence']:
            sequence = info['sequence']
            gaps_flag = 'no'
        else:
            continue
        if sequence not in seen_sequences:
            seen_sequences.add(sequence)
            tag = f"{pdb_id}_{chain_id}"
            header_parts = [f"{tag}", f"entity_id={info['entity_id']}", f"gaps={gaps_flag}"]
            db_refs = info.get('db_references', [])
            uniprot_refs = [ref for ref in db_refs if ref['db_name'] == 'UniProt']
            other_refs = [ref for ref in db_refs if ref['db_name'] != 'UniProt']
            for ref in uniprot_refs + other_refs:
                header_parts.append(f"db={ref['db_name']}")
                header_parts.append(f"id={ref['db_code']}")
            header = " ".join(header_parts)
            with open(outdir / f"{tag}_seqs.fasta", "w") as f:
                f.write(f">{header}\n")
                f.write(sequence + "\n")
            written_count += 1
    return written_count

def write_atom_coordinates_pdb(chains_info, pdb_id, outdir, include_hetatm=False, include_nonprotein=False):
    for chain_id, info in chains_info.items():
        if info.get('entity_type') != 'polymer' and not include_nonprotein:
            continue
        tag = f"{pdb_id}_{chain_id}"
        atoms_path = outdir / f"{tag}_atoms.pdb"
        with open(atoms_path, "w") as f:
            for atom_num, atom_data in enumerate(info.get('all_atoms', []), start=1):
                hetflag = atom_data.get('hetflag', '').strip()
                is_het = hetflag not in ('', ' ')
                if is_het and not include_hetatm:
                    continue
                record_type = 'HETATM' if is_het else 'ATOM'
                raw_name = atom_data['atom_name']
                element = atom_data['element']
                if len(raw_name) < 4 and len(element) == 1:
                    atom_name_field = f" {raw_name:<3s}"
                else:
                    atom_name_field = f"{raw_name:<4s}"
                resname_field = f"{atom_data['resname']:>3s}"
                seq_num = atom_data['seq_num']
                icode = atom_data['icode'] if atom_data['icode'] else ' '
                chain_field = chain_id[0] if chain_id else ' '
                f.write(
                    f"{record_type:<6s}"
                    f"{atom_num:>5d} "
                    f"{atom_name_field}"
                    f" "
                    f"{resname_field} "
                    f"{chain_field}"
                    f"{seq_num:>4d}"
                    f"{icode}"
                    f"   "
                    f"{atom_data['x']:>8.3f}"
                    f"{atom_data['y']:>8.3f}"
                    f"{atom_data['z']:>8.3f}"
                    f"{atom_data['occupancy']:>6.2f}"
                    f"{atom_data['b_factor']:>6.2f}"
                    f"          "
                    f"{element:>2s}"
                    f"\n"
                )
            f.write("END\n")

def write_atom_coordinates_tsv(chains_info, pdb_id, outdir):
    for chain_id, info in chains_info.items():
        tag = f"{pdb_id}_{chain_id}"
        atoms_path = outdir / f"{tag}_atoms.tsv"
        with open(atoms_path, "w") as f:
            f.write("record\tatom_num\tatom_name\tresname\tchain\tres_seq\t"
                    "x\ty\tz\toccupancy\tb_factor\telement\n")
            for atom_num, atom_data in enumerate(info.get('all_atoms', []), start=1):
                hetflag = atom_data.get('hetflag', '').strip()
                record = 'HETATM' if hetflag not in ('', ' ') else 'ATOM'
                f.write(f"{record}\t{atom_num}\t{atom_data['atom_name']}\t"
                        f"{atom_data['resname']}\t{chain_id}\t{atom_data['auth_seq']}\t"
                        f"{atom_data['x']:.3f}\t{atom_data['y']:.3f}\t{atom_data['z']:.3f}\t"
                        f"{atom_data['occupancy']:.2f}\t{atom_data['b_factor']:.2f}\t"
                        f"{atom_data['element']}\n")

def main():
    parser = argparse.ArgumentParser(description='Download mmCIF and extract per-chain PDB/TSV files')
    parser.add_argument('pdb_id', help='PDB ID to download')
    parser.add_argument('output_directory', help='Output directory')
    parser.add_argument('--format', choices=['tsv', 'pdb'], default='tsv',
                        help='Output format for atom coordinates (default: tsv)')
    parser.add_argument('--include-hetatm', action='store_true', default=False,
                        help='Include HETATM records (default: False)')
    parser.add_argument('--include-nonprotein', action='store_true', default=False,
                        help='Include non-protein chains (default: False)')
    args = parser.parse_args()
    pdb_id = args.pdb_id.strip().upper()
    outdir = Path(args.output_directory).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    mmcif_path = download_mmcif(pdb_id, outdir)
    chains = parse_with_biopython(mmcif_path, pdb_id)
    if not chains:
        print("No chains found.")
        return
    write_fasta_files(chains, pdb_id, outdir)
    if args.format == 'tsv':
        write_atom_coordinates_tsv(chains, pdb_id, outdir)
    else:
        write_atom_coordinates_pdb(chains, pdb_id, outdir, args.include_hetatm, args.include_nonprotein)

if __name__ == "__main__":
    main()


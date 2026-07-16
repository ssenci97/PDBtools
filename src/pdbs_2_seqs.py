#!/usr/bin/env python3
import argparse
import os
import pandas as pd
from pathlib import Path
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1


def get_pdb_sequence(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    sequence_parts = []

    if len(structure) > 0:
        first_model = structure[0]

        for chain in first_model:
            for residue in chain:
                if residue.id[0] == " ":
                    try:
                        sequence_parts.append(seq1(residue.get_resname()))
                    except Exception:
                        continue

    return "".join(sequence_parts)


def main():
    parser = argparse.ArgumentParser(description="Extract sequences from paired PDB files in subdirectories and save to CSV/FASTA.")
    parser.add_argument("-i", "--input", required=True, type=str, help="Path to the base input directory containing subdirectories with PDB files.")
    parser.add_argument("-o", "--output", required=True, type=str, help="Path to the output CSV file.")
    parser.add_argument("--tofasta", action="store_true", help="Also save individual FASTA files alongside the CSV.")
    args = parser.parse_args()

    base_dir = Path(args.input)
    output_path = Path(args.output)

    if not base_dir.is_dir():
        print(f"Error: Input directory '{base_dir}' does not exist.")
        return

    pdb_files = list(base_dir.rglob("*.pdb"))

    if not pdb_files:
        print(f"No .pdb files found in {base_dir} or its subdirectories.")
        return

    print(f"Found {len(pdb_files)} PDB files. Processing...")
    data = []

    for pdb_file in pdb_files:
        fname = pdb_file.stem
        complex_name = pdb_file.parent.name

        try:
            sequence = get_pdb_sequence(str(pdb_file))
            parts = fname.split("_", 1)
            pdb_id = parts[0] if len(parts) > 0 else fname
            chain = parts[1] if len(parts) > 1 else ""
            entry_id = f"{pdb_id}_{chain}" if chain else pdb_id

            data.append({
                "complex": complex_name,
                "id": entry_id,
                "sequence": sequence
            })
        except Exception as e:
            print(f"Error processing {pdb_file}: {e}")

    output_dir = output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame(data, columns=["complex", "id", "sequence"]).drop_duplicates()
    df.to_csv(output_path, index=False)
    print(f"Successfully saved dataframe to {output_path}")

    if args.tofasta:
        for entry in data:
            fasta_name = f"{entry['id']}.fasta"
            fasta_path = output_dir / fasta_name

            with open(fasta_path, "w") as f:
                f.write(f">{entry['id']}\n")
                f.write(f"{entry['sequence']}\n")
        print(f"Saved {len(data)} FASTA files to {output_dir}")
    else:
        print("Warning: FASTA files were not saved. Add --tofasta to also generate individual FASTA files.")


if __name__ == "__main__":
    main()


#!/usr/bin/env python3
"""Align PDB structures pairwise between two file lists.

Usage:
    python align_pdbs.py --f1 list1.txt --f2 list2.txt --output results.tsv

Each line in --f1 and --f2 must contain one PDB file path.
Every file in --f1 is paired with every file in --f2 (unordered pairs).
"""

import os
import sys
import argparse
import hashlib
import gzip
import shutil
import tempfile
import pandas as pd
from pymol import cmd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Pairwise RMSD alignment of PDB files."
    )
    parser.add_argument("--f1", required=True, help="File with one PDB path per line")
    parser.add_argument("--f2", required=True, help="File with one PDB path per line")
    parser.add_argument("--output", required=True, help="Output TSV file path")
    return parser.parse_args()


def load_file_list(path: str) -> list[str]:
    """Load a list of file paths from a text file (one per line)."""
    with open(path, "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    return [os.path.abspath(line) for line in lines]


def align_pdb_structures(mol1_path: str, mol2_path: str) -> pd.DataFrame:
    """Align two PDB structures with PyMOL and return metrics as a DataFrame."""
    with tempfile.TemporaryDirectory() as tmp_dir:

        def prepare(path: str) -> str:
            if path.endswith(".gz"):
                h = hashlib.md5(path.encode()).hexdigest()[:8]
                out = os.path.join(tmp_dir, f"{h}_{os.path.basename(path)[:-3]}")
                with gzip.open(path, "rb") as f_in, open(out, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
                return out
            return path

        mol1_local, mol2_local = prepare(mol1_path), prepare(mol2_path)
        mol1_name = f"mol1_{hashlib.md5(mol1_path.encode()).hexdigest()[:8]}"
        mol2_name = f"mol2_{hashlib.md5(mol2_path.encode()).hexdigest()[:8]}"

        cmd.reinitialize()
        cmd.load(mol1_local, mol1_name)
        cmd.load(mol2_local, mol2_name)

        (rmsd_after, n_atoms_aligned, n_cycles, rmsd_before,
         n_atoms_pre, score, n_res_aligned) = cmd.align(mol1_name, mol2_name)

        df = pd.DataFrame([{
            "mol_1": mol1_path, "mol_2": mol2_path,
            "rmsd": rmsd_after, "rmsd_before_refinement": rmsd_before,
            "n_atoms_aligned": n_atoms_aligned, "n_cycles": n_cycles,
            "n_atoms_pre_refinement": n_atoms_pre, "score": score,
            "n_residues_aligned": n_res_aligned,
        }])

    return df


def main():
    args = parse_args()

    outpath = os.path.abspath(args.output)

    if os.path.exists(outpath):
        print(f"[SKIP] Output file already exists: {outpath}")
        sys.exit(0)

    files1 = load_file_list(args.f1)
    files2 = load_file_list(args.f2)

    if not files1:
        print(f"[ERROR] No file paths loaded from {args.f1}")
        sys.exit(1)
    if not files2:
        print(f"[ERROR] No file paths loaded from {args.f2}")
        sys.exit(1)

    # Sanity-check: warn about missing files
    for f in files1:
        if not os.path.exists(f):
            print(f"[WARN] File not found (f1): {f}")
    for f in files2:
        if not os.path.exists(f):
            print(f"[WARN] File not found (f2): {f}")

    # Every file in f1 paired with every file in f2.
    # mol_1 always from f1, mol_2 always from f2.
    rows: list[pd.DataFrame] = []
    for p1 in files1:
        for p2 in files2:
            print(f"  {os.path.basename(p1)}  <->  {os.path.basename(p2)}")
            df_row = align_pdb_structures(mol1_path=p1, mol2_path=p2)
            rows.append(df_row)

    if not rows:
        print("[WARN] No pairs to align.")
        sys.exit(0)

    os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
    final_df = pd.concat(rows, ignore_index=True)
    final_df.to_csv(outpath, sep="\t", index=False)
    print(f"[DONE] Results written to {outpath}  ({len(final_df)} rows)")


if __name__ == "__main__":
    main()

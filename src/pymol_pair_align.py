#!/usr/bin/env python3
import os
import gzip
import shutil
import tempfile
import hashlib
import argparse
import pandas as pd
import pymol
from pymol import cmd
import itertools


def _prepare_pdb_path(pdb_path: str, tmp_dir: str) -> str:
    if pdb_path.endswith(".gz"):
        base_name = os.path.basename(pdb_path)[:-3]
        path_hash = hashlib.md5(pdb_path.encode()).hexdigest()[:8]
        decompressed_path = os.path.join(tmp_dir, f"{path_hash}_{base_name}")
        with gzip.open(pdb_path, 'rb') as f_in:
            with open(decompressed_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        return decompressed_path
    return pdb_path


def _unique_obj_name(pdb_path: str, prefix: str) -> str:
    h = hashlib.md5(pdb_path.encode()).hexdigest()[:8]
    return f"{prefix}_{h}"


def align_pdb_structures(mol1_path: str, mol2_path: str) -> pd.DataFrame:
    with tempfile.TemporaryDirectory() as tmp_dir:
        mol1_local_path = _prepare_pdb_path(mol1_path, tmp_dir)
        mol2_local_path = _prepare_pdb_path(mol2_path, tmp_dir)

        mol1_name = _unique_obj_name(mol1_path, "mol1")
        mol2_name = _unique_obj_name(mol2_path, "mol2")

        cmd.reinitialize()
        cmd.load(mol1_local_path, mol1_name)
        cmd.load(mol2_local_path, mol2_name)

        result = cmd.align(mol1_name, mol2_name)

        (rmsd_after_refinement,
         n_atoms_aligned,
         n_cycles,
         rmsd_before_refinement,
         n_atoms_pre_refinement,
         score,
         n_residues_aligned) = result

        df = pd.DataFrame([{
            "mol_1": mol1_path,
            "mol_2": mol2_path,
            "rmsd": rmsd_after_refinement,
            "rmsd_before_refinement": rmsd_before_refinement,
            "n_atoms_aligned": n_atoms_aligned,
            "n_cycles": n_cycles,
            "n_atoms_pre_refinement": n_atoms_pre_refinement,
            "score": score,
            "n_residues_aligned": n_residues_aligned
        }])

    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required=True)
    parser.add_argument("--outpath", required=True)
    args = parser.parse_args()

    pymol.pymol_argv = ['pymol', '-qc']
    pymol.finish_launching(['pymol', '-qc'])

    out_dir = os.path.dirname(args.outpath)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    pdb_files = []

    for root, dirs, files in os.walk(args.input_dir):
        if os.path.basename(root) == '3_emref':
            for f in files:
                if f.endswith('.pdb.gz'):
                    pdb_files.append(os.path.join(root, f))

    all_results = []

    for path_1, path_2 in itertools.combinations(pdb_files, 2):
        df_result = align_pdb_structures(path_1, path_2)
        all_results.append(df_result)
        print("Done one\n")

    combined_df = pd.concat(all_results, ignore_index=True)
    combined_df.to_csv(args.outpath, sep='\t', index=False)

    cmd.quit()


if __name__ == "__main__":
    main()
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


#### FUNCTIONS
def align_pdb_structures(mol1_path: str, mol2_path: str, return_structure: bool = False):
    with tempfile.TemporaryDirectory() as tmp_dir:

        def prepare(path):
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

        if return_structure:
            aligned_structures = {
                "mol_1": cmd.get_pdbstr(mol1_name),
                "mol_2": cmd.get_pdbstr(mol2_name),
            }
            return df, aligned_structures

    return df

##########################
def plot_aligned_structures(aligned_structures, width=800, height=600):
    import py3Dmol

    def get_chains(pdb_str):
        chains = []
        for line in pdb_str.splitlines():
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21]
                if chain not in chains:
                    chains.append(chain)
        return chains

    hot_colors = ["orange", "yellow", "red", "gold"]
    cold_colors = ["lightblue", "green", "teal", "cyan"]

    view = py3Dmol.view(width=width, height=height)

    view.addModel(aligned_structures["mol_1"], "pdb")
    for i, chain in enumerate(get_chains(aligned_structures["mol_1"])):
        view.setStyle({"model": 0, "chain": chain}, {"cartoon": {"color": hot_colors[i % len(hot_colors)]}})

    view.addModel(aligned_structures["mol_2"], "pdb")
    for i, chain in enumerate(get_chains(aligned_structures["mol_2"])):
        view.setStyle({"model": 1, "chain": chain}, {"cartoon": {"color": cold_colors[i % len(cold_colors)]}})

    view.zoomTo()
    return view

#######################################
def plot_rmsd_heatmap(aln_df, rmsd_col='rmsd', save_path=None, label_replace=None, title=None):
    pdb_ids = sorted(set(aln_df['pdb_id1']) | set(aln_df['pdb_id2']))
    n = len(pdb_ids)

    matrix = np.full((n, n), np.nan)
    id_to_idx = {pdb_id: i for i, pdb_id in enumerate(pdb_ids)}

    for _, row in aln_df.iterrows():
        i = id_to_idx[row['pdb_id1']]
        j = id_to_idx[row['pdb_id2']]
        matrix[i, j] = row[rmsd_col]
        matrix[j, i] = row[rmsd_col]

    mask = np.triu(np.ones_like(matrix, dtype=bool), k=1)

    fig, ax = plt.subplots(figsize=(max(8, n * 0.6), max(6, n * 0.5)))

    cmap = mcolors.LinearSegmentedColormap.from_list(
        'rmsd_cmap', ['#006400', '#228B22', '#9ACD32', '#D3D3D3', '#696969']
    )

    sns.heatmap(
        matrix, mask=mask, xticklabels=pdb_ids, yticklabels=pdb_ids,
        cmap=cmap, annot=True, fmt='.1f', linewidths=0.5, square=True,
        cbar_kws={'label': 'RMSD (Å)'}, ax=ax, vmin=0, vmax=5
    )

    if label_replace:
        xlabels = [l.get_text() for l in ax.get_xticklabels()]
        ylabels = [l.get_text() for l in ax.get_yticklabels()]
        for rep in label_replace:
            xlabels = [x.replace(rep, '') for x in xlabels]
            ylabels = [y.replace(rep, '') for y in ylabels]
        ax.set_xticklabels(xlabels)
        ax.set_yticklabels(ylabels)

    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title(title or 'RMSD All-vs-All Triangular Heatmap')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()

    return fig, ax

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

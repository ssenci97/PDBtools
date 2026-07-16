#!/usr/bin/env python3
import mdtraj as md
import numpy as np
import pandas as pd
import argparse
import os
import sys

# --- TUNABLES ---
DEFAULT_CUTOFF_NM = 0.8
DEFAULT_MODE = "inter"
DEFAULT_OUTPUT = "contacts.tsv"


def main():
    parser = argparse.ArgumentParser(description="Extract residue contacts from an mmCIF file.")
    parser.add_argument("input", help="Path to the input mmCIF file (.cif)")
    parser.add_argument("-o", "--output", default=DEFAULT_OUTPUT, help="Output TSV filename")
    parser.add_argument("-c", "--cutoff", type=float, default=DEFAULT_CUTOFF_NM, help="Distance cutoff in nanometers")
    parser.add_argument("--mode", choices=["intra", "inter", "both"], default=DEFAULT_MODE, help="Contact mode")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: File {args.input} not found.")
        sys.exit(1)

    print(f"Loading {args.input}...")
    try:
        traj = md.load(args.input)
    except Exception as e:
        print(f"Error loading file: {e}")
        sys.exit(1)

    top = traj.topology
    pdb_id = os.path.splitext(os.path.basename(args.input))[0]

    ca_atoms = [a for a in top.atoms if a.name == 'CA']
    if len(ca_atoms) < 2:
        print("Error: Fewer than 2 CA atoms found.")
        sys.exit(1)

    pairs = []
    for i, ai in enumerate(ca_atoms):
        for j, aj in enumerate(ca_atoms[i + 1:], start=i + 1):
            same_chain = ai.residue.chain.index == aj.residue.chain.index
            if args.mode == "intra" and not same_chain:
                continue
            if args.mode == "inter" and same_chain:
                continue
            pairs.append((ai.index, aj.index))

    pairs = np.array(pairs)
    print(f"Computing {len(pairs)} CA-CA distances (mode={args.mode}, cutoff={args.cutoff} nm)...")
    if len(pairs) == 0:
        print("No pairs found.")
        sys.exit(0)

    dists = md.compute_distances(traj, pairs)[0]
    mask = dists <= args.cutoff
    pairs = pairs[mask]
    dists = dists[mask]
    print(f"Kept {len(pairs)} pairs within cutoff")

    contact_data = []
    for (ai_idx, aj_idx), d in zip(pairs, dists):
        ri = top.atom(ai_idx).residue
        rj = top.atom(aj_idx).residue
        contact_data.append({
            'pdb_id': pdb_id,
            'res_i': ri.index,
            'res_j': rj.index,
            'chain_i': ri.chain.chain_id,
            'chain_j': rj.chain.chain_id,
            'resname_i': ri.name,
            'resname_j': rj.name,
            'resseq_i': ri.resSeq,
            'resseq_j': rj.resSeq,
            'aminoacid_i': f"{ri.name}{ri.resSeq}",
            'aminoacid_j': f"{rj.name}{rj.resSeq}",
            'dist_nm': round(d, 6),
            'dist_A': round(d * 10, 3),
        })

    df = pd.DataFrame(contact_data)
    df.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote {len(pairs)} contacts to {args.output}")


if __name__ == "__main__":
    main()


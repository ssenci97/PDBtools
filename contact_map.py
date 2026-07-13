#!/usr/bin/env python3
import argparse
import mdtraj as md
import numpy as np

def main():
    parser = argparse.ArgumentParser(description='Calculate inter-chain CA contact map from PDB')
    parser.add_argument('pdb_file', help='Input PDB file')
    parser.add_argument('-o', '--output', default='contacts.tsv', help='Output TSV file')
    parser.add_argument('--cutoff', type=float, default=None, help='Contact cutoff in nm')

    args = parser.parse_args()

    traj = md.load(args.pdb_file)
    top = traj.topology

    ca_atoms = [a for a in top.atoms if a.name == 'CA']
    print(f"Found {len(ca_atoms)} CA atoms")

    if len(ca_atoms) < 2:
        print("Error: Fewer than 2 CA atoms found.")
        return

    pairs = []
    for i, ai in enumerate(ca_atoms):
        for j, aj in enumerate(ca_atoms[i+1:], start=i+1):
            if ai.residue.chain.index != aj.residue.chain.index:
                pairs.append((ai.index, aj.index))

    pairs = np.array(pairs)
    print(f"Computing {len(pairs)} inter-chain CA-CA distances...")

    if len(pairs) == 0:
        print("No inter-chain pairs found.")
        return

    dists = md.compute_distances(traj, pairs)[0]

    if args.cutoff is not None:
        mask = dists <= args.cutoff
        pairs = pairs[mask]
        dists = dists[mask]
        print(f"Kept {len(pairs)} pairs within {args.cutoff} nm cutoff")

    with open(args.output, 'w') as f:
        f.write("res_i\tres_j\tchain_i\tchain_j\tresname_i\tresname_j\t"
                "resseq_i\tresseq_j\tdist_nm\tdist_A\n")

        for (ai_idx, aj_idx), d in zip(pairs, dists):
            ri = top.atom(ai_idx).residue
            rj = top.atom(aj_idx).residue

            f.write(f"{ri.index}\t{rj.index}\t"
                    f"{ri.chain.chain_id}\t{rj.chain.chain_id}\t"
                    f"{ri.name}\t{rj.name}\t"
                    f"{ri.resSeq}\t{rj.resSeq}\t"
                    f"{d:.6f}\t{d*10:.6f}\n")

    print(f"Wrote {len(pairs)} inter-chain contacts to {args.output}")

if __name__ == '__main__':
    main()
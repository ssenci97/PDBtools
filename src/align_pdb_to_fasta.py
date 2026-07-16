#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
from Bio.PDB import PDBParser
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import pandas as pd
# =============================================================================
# TUNABLES
# =============================================================================
ALIGNMENT_MATCH_SCORE = 1.0
ALIGNMENT_MISMATCH_SCORE = -1.0
ALIGNMENT_OPEN_GAP_SCORE = -2.0
ALIGNMENT_EXTEND_GAP_SCORE = -0.5
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
                        resname = residue.get_resname().strip().upper()
                        if resname in protein_letters_3to1_extended:
                            sequence_parts.append(protein_letters_3to1_extended[resname])
                    except Exception:
                        continue
    return "".join(sequence_parts)
def load_fasta_sequence(fasta_path):
    return str(SeqIO.read(fasta_path, "fasta").seq)
def run_alignment(ref_seq, query_seq, mode="global"):
    aligner = PairwiseAligner(
        mode=mode,
        match_score=ALIGNMENT_MATCH_SCORE,
        mismatch_score=ALIGNMENT_MISMATCH_SCORE,
        open_gap_score=ALIGNMENT_OPEN_GAP_SCORE,
        extend_gap_score=ALIGNMENT_EXTEND_GAP_SCORE,
    )
    alignments = aligner.align(ref_seq, query_seq)
    if not alignments:
        return None
    return alignments[0]
def get_alignment_start_offsets(alignment):
    """
    Return the 0-based starting index of the aligned block in the ORIGINAL
    ref/query sequences. In 'global' mode this is normally 0, but in 'local'
    mode Biopython's alignment[0]/alignment[1] strings only cover the aligned
    span (not padded to the full sequence), so without this offset all
    downstream numbering silently resets to 1 instead of reflecting the true
    position in the full-length sequence.
    """
    try:
        ref_blocks, query_blocks = alignment.aligned
        if len(ref_blocks) > 0:
            return int(ref_blocks[0][0]), int(query_blocks[0][0])
    except Exception:
        pass
    return 0, 0
def analyze_fragments_and_range(a, b, ref_start_offset=0, query_start_offset=0):
    """
    Walk the aligned pair (a = ref/UniProt aligned string, b = query/PDB aligned string)
    and derive per-fragment statistics. A "fragment" is a contiguous stretch of the
    query (PDB) sequence with no gaps in the query. For each fragment we record both
    the UniProt (ref) numbering range and the PDB (query) numbering range separately,
    since gaps on either side can make these ranges diverge in length/position.
    ref_start_offset / query_start_offset seed the position counters at the true
    0-based start of the aligned block within the full-length sequences (needed
    because local-mode alignment strings do not start at position 0 of the
    original sequence).
    """
    frag_lengths = []
    ref_ranges = []
    query_ranges = []
    gap_lengths = []
    ref_pos = ref_start_offset
    query_pos = query_start_offset
    current_query_start = None
    current_ref_start = None
    current_length = 0
    current_gap = 0
    for x, y in zip(a, b):
        if x != "-":
            ref_pos += 1
        if y != "-":
            query_pos += 1
        if y != "-":
            if current_query_start is None:
                if current_gap > 0:
                    gap_lengths.append(current_gap)
                current_query_start = query_pos
                current_ref_start = ref_pos
                current_length = 0
                current_gap = 0
            current_length += 1
        else:
            if current_query_start is not None:
                frag_lengths.append(current_length)
                ref_end = ref_pos - 1 if x == "-" else ref_pos
                query_end = query_pos
                ref_ranges.append(f"{current_ref_start}-{ref_end}")
                query_ranges.append(f"{current_query_start}-{query_end}")
                current_query_start = None
                current_length = 0
                current_gap = 1
            else:
                current_gap += 1
    if current_query_start is not None:
        frag_lengths.append(current_length)
        ref_end = ref_pos
        query_end = query_pos
        ref_ranges.append(f"{current_ref_start}-{ref_end}")
        query_ranges.append(f"{current_query_start}-{query_end}")
    if current_gap > 0:
        gap_lengths.append(current_gap)
    frag_lengths_str = ",".join(map(str, frag_lengths)) if frag_lengths else ""
    frag_ranges_upkb_str = ";".join(ref_ranges) if ref_ranges else ""
    frag_ranges_pdb_str = ";".join(query_ranges) if query_ranges else ""
    gap_lengths_str = ",".join(map(str, gap_lengths)) if gap_lengths else ""
    ref_range = ""
    if ref_ranges:
        starts = [int(r.split("-")[0]) for r in ref_ranges]
        ends = [int(r.split("-")[1]) for r in ref_ranges]
        ref_range = f"{min(starts)}-{max(ends)}"
    return {
        "num_frags": len(frag_lengths),
        "frag_lengths": frag_lengths_str,
        "frag_ranges_upkb": frag_ranges_upkb_str,
        "frag_ranges_pdb": frag_ranges_pdb_str,
        "gap_lengths": gap_lengths_str,
        "ref_range": ref_range,
    }
def compute_alignment_stats(ref_seq, query_seq, alignment, mode):
    a = alignment[0]
    b = alignment[1]
    ref_start_offset, query_start_offset = get_alignment_start_offsets(alignment)
    matches = mismatches = insertions = deletions = 0
    for x, y in zip(a, b):
        if x == "-" and y == "-":
            continue
        elif x == "-":
            insertions += 1
        elif y == "-":
            deletions += 1
        elif x == y:
            matches += 1
        else:
            mismatches += 1
    aligned_length = matches + mismatches + insertions + deletions
    query_aligned_len = len(query_seq) - deletions
    frag_stats = analyze_fragments_and_range(a, b, ref_start_offset, query_start_offset)
    identity_percent = round(matches / aligned_length * 100, 2) if aligned_length > 0 else 0.0
    coverage_percent = round(query_aligned_len / len(query_seq) * 100, 2) if len(query_seq) > 0 else 0.0
    score = round(
        matches
        - 0.5 * mismatches
        - 1.5 * (insertions + deletions),
        2,
    )
    stats = {
        "ref_length": len(ref_seq),
        "query_length": len(query_seq),
        "aligned_length": aligned_length,
        "matches": matches,
        "mismatches": mismatches,
        "insertions_in_query": insertions,
        "deletions_in_query": deletions,
        **frag_stats,
        "identity_percent": identity_percent,
        "query_coverage_percent": coverage_percent,
        "alignment_score": score,
        "alignment_mode": mode,
    }
    return stats, a, b, ref_start_offset, query_start_offset
def build_ranges_summary(stats, ref_id="ref", query_id="query"):
    """Build a single-row DataFrame with all alignment statistics."""
    row = {
        "ref_id": ref_id,
        "query_id": query_id,
        **stats,
    }
    return pd.DataFrame([row])
def build_residues_dataframe(ref_aligned, query_aligned, ref_id="ref", query_id="query",
                              ref_start_offset=0, query_start_offset=0):
    rows = []
    ref_pos = ref_start_offset
    query_pos = query_start_offset
    alignment_pos = 0
    for x, y in zip(ref_aligned, query_aligned):
        alignment_pos += 1
        if x != "-":
            ref_pos += 1
        if y != "-":
            query_pos += 1
        if x == "-" and y == "-":
            continue
        rows.append(
            {
                "alignment_position": alignment_pos,
                "ref_id": ref_id,
                "ref_position": ref_pos if x != "-" else -1,
                "ref_aa": x,
                "query_id": query_id,
                "query_position": query_pos if y != "-" else -1,
                "query_aa": y,
                "match": "1" if (x == y and x != "-") else "0",
            }
        )
    return pd.DataFrame(rows)
def debug_pdb_structure(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    print(f"\n=== PDB STRUCTURE DEBUG ===")
    print(f"Number of models: {len(structure)}")
    if len(structure) == 0:
        print("No models found!")
        return
    model = structure[0]
    print(f"Number of chains: {len(list(model.get_chains()))}")
    for i, chain in enumerate(model):
        print(f"\nChain {i+1}: ID='{chain.id}'")
        residues = list(chain.get_residues())
        print(f"  Number of residues: {len(residues)}")
        for j, residue in enumerate(residues[:20]):
            resname = residue.get_resname().strip().upper()
            het_flag = residue.id[0]
            resnum = residue.id[1]
            icode = residue.id[2]
            print(f"    Residue {j+1}: {resname} {resnum}{icode} (het_flag='{het_flag}')")
        if len(residues) > 20:
            print(f"    ... and {len(residues) - 20} more residues")
def main():
    parser = argparse.ArgumentParser(
        description="Align PDB-derived sequence to a FASTA reference sequence."
    )
    parser.add_argument("--PDB", required=True, help="Path to PDB structure file")
    parser.add_argument("--FASTA", required=True, help="Path to FASTA reference sequence")
    parser.add_argument(
        "--mode",
        required=True,
        choices=["local", "global"],
        help="Alignment mode: local or global",
    )
    parser.add_argument(
        "--output-ranges",
        default=None,
        help="Output TSV path for alignment summary statistics",
    )
    parser.add_argument(
        "--output-residues",
        default=None,
        help="Output TSV path for per-residue alignment details",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print debug information about PDB structure",
    )
    args = parser.parse_args()
    pdb_path = Path(args.PDB)
    fasta_path = Path(args.FASTA)
    if not pdb_path.exists():
        print(f"Error: PDB file not found: {pdb_path}", file=sys.stderr)
        sys.exit(1)
    if not fasta_path.exists():
        print(f"Error: FASTA file not found: {fasta_path}", file=sys.stderr)
        sys.exit(1)
    if args.debug:
        debug_pdb_structure(str(pdb_path))
    print(f"\nLoading PDB sequence from {pdb_path} ...")
    pdb_seq = get_pdb_sequence(str(pdb_path))
    print(f"PDB sequence length: {len(pdb_seq)}")
    if len(pdb_seq) > 0:
        print(f"First 20 residues: {pdb_seq[:20]}")
    print(f"\nLoading FASTA sequence from {fasta_path} ...")
    fasta_seq = load_fasta_sequence(str(fasta_path))
    print(f"FASTA sequence length: {len(fasta_seq)}")
    if not pdb_seq:
        print("Error: No sequence could be extracted from the PDB file.", file=sys.stderr)
        sys.exit(1)
    if not fasta_seq:
        print("Error: No sequence could be read from the FASTA file.", file=sys.stderr)
        sys.exit(1)
    print(f"\nRunning {args.mode} alignment ...")
    alignment = run_alignment(fasta_seq, pdb_seq, mode=args.mode)
    if alignment is None:
        print("Error: Alignment failed.", file=sys.stderr)
        sys.exit(1)
    stats, ref_aligned, query_aligned, ref_start_offset, query_start_offset = compute_alignment_stats(
        fasta_seq, pdb_seq, alignment, args.mode
    )
    print("\n" + "=" * 60)
    print("ALIGNMENT SUMMARY")
    print("=" * 60)
    for key, value in stats.items():
        print(f"  {key:<25}: {value}")
    print("=" * 60)
    if args.output_ranges:
        ranges_df = build_ranges_summary(
            stats, ref_id=fasta_path.stem, query_id=pdb_path.stem
        )
        ranges_path = Path(args.output_ranges)
        ranges_path.parent.mkdir(parents=True, exist_ok=True)
        ranges_df.to_csv(ranges_path, sep="\t", index=False)
        print(f"\nSaved alignment summary -> {ranges_path}")
        print(ranges_df.to_string(index=False))
    if args.output_residues:
        residues_df = build_residues_dataframe(
            ref_aligned, query_aligned, ref_id=fasta_path.stem, query_id=pdb_path.stem,
            ref_start_offset=ref_start_offset, query_start_offset=query_start_offset,
        )
        residues_path = Path(args.output_residues)
        residues_path.parent.mkdir(parents=True, exist_ok=True)
        residues_df.to_csv(residues_path, sep="\t", index=False)
        print(f"\nSaved per-residue alignment -> {residues_path}")
        print(f"Total aligned positions: {len(residues_df)}")
if __name__ == "__main__":
    main()

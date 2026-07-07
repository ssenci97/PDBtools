#!/usr/bin/env python3
from pathlib import Path
from Bio.PDB import PDBParser, PPBuilder
from Bio.Align import PairwiseAligner
import pandas as pd


def get_pdb_sequence(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('X', pdb_file)
    return ''.join(str(pp.get_sequence()) for pp in PPBuilder().build_peptides(structure))


def analyze_fragments_and_range(a, b):
    """
    Query fragments (continuous in PDB), but ranges projected on REFERENCE.
    gap_lengths = gaps between these fragments (in reference coordinates).
    """
    query_fragments = []
    ref_ranges = []           # projected on reference
    gap_lengths = []

    ref_pos = 0
    query_pos = 0

    current_query_start = None
    current_ref_start = None
    current_length = 0
    current_gap = 0

    for x, y in zip(a, b):          # a = ref, b = query
        if x != '-': ref_pos += 1
        if y != '-': query_pos += 1

        if y != '-':                                 # Query residue (PDB)
            if current_query_start is None:
                # Start new query fragment
                if current_gap > 0:
                    gap_lengths.append(current_gap)
                current_query_start = query_pos
                current_ref_start = ref_pos
                current_length = 0
                current_gap = 0

            current_length += 1
        else:                                        # Gap in query
            if current_query_start is not None:
                # End current fragment
                query_fragments.append(current_length)
                # Range on reference
                ref_end = ref_pos - 1 if x == '-' else ref_pos
                ref_ranges.append(f"{current_ref_start}-{ref_end}")
                current_query_start = None
                current_length = 0
                current_gap = 1
            else:
                current_gap += 1

    # Last fragment
    if current_query_start is not None:
        query_fragments.append(current_length)
        ref_end = ref_pos
        ref_ranges.append(f"{current_ref_start}-{ref_end}")

    # Last gap
    if current_gap > 0:
        gap_lengths.append(current_gap)

    frag_lengths_str = ','.join(map(str, query_fragments)) if query_fragments else ''
    frag_ranges_str = ';'.join(ref_ranges) if ref_ranges else ''
    gap_lengths_str = ','.join(map(str, gap_lengths)) if gap_lengths else ''

    ref_range = f"{min([int(r.split('-')[0]) for r in ref_ranges])}-{max([int(r.split('-')[1]) for r in ref_ranges])}" \
                if ref_ranges else ''

    return {
        'num_fragments': len(query_fragments),
        'fragment_lengths': frag_lengths_str,      # length in query
        'fragment_ranges': frag_ranges_str,        # on reference!
        'gap_lengths': gap_lengths_str,
        'ref_range': ref_range
    }


def score_alignment(ref_seq, query_seq):
    aligner = PairwiseAligner(
        mode='local',
        match_score=1,
        mismatch_score=-1,
        open_gap_score=-2,
        extend_gap_score=-0.5
    )
   
    alignments = aligner.align(ref_seq, query_seq)
    if not alignments:
        return None
   
    a, b = alignments[0]

    matches = mismatches = insertions = deletions = 0
    for x, y in zip(a, b):
        if x == '-' and y == '-': continue
        elif x == '-': insertions += 1
        elif y == '-': deletions += 1
        elif x == y: matches += 1
        else: mismatches += 1

    aligned_length = matches + mismatches + insertions + deletions
    query_aligned_len = len(query_seq) - deletions

    extra = analyze_fragments_and_range(a, b)

    return {
        'file': None,
        'query_length': len(query_seq),
        'aligned_length': aligned_length,
        'matches': matches,
        'mismatches': mismatches,
        'insertions_in_query': insertions,
        'deletions_in_query': deletions,
        **extra,
        'identity_percent': round(matches / aligned_length * 100, 2) if aligned_length > 0 else 0,
        'query_coverage_percent': round(query_aligned_len / len(query_seq) * 100, 2),
        'alignment_score': round(matches - 0.5 * mismatches - 1.5 * (insertions + deletions), 2),
    }


def main(fasta_uniprot, fasta_dir, output):
    ref_text = Path(fasta_uniprot).read_text().strip()
    ref_seq = ''.join(ref_text.split('\n')[1:]).replace('\n', '')
    stem = Path(fasta_uniprot).stem
    rows = []

    for pdb in sorted(Path(fasta_dir).glob(f'{stem}*.pdb')):
        query_seq = get_pdb_sequence(pdb)
        stats = score_alignment(ref_seq, query_seq)
        if stats:
            stats['file'] = pdb.name
            rows.append(stats)
           
            print(f"{pdb.name}: {stats['identity_percent']}% id | "
                  f"{stats['num_fragments']} frags {stats['fragment_lengths']} | "
                  f"gaps:{stats['gap_lengths']} | ranges(ref): {stats['fragment_ranges']}")
        else:
            print(f"{pdb.name}: No alignment")

    if not rows:
        print("No alignments found.")
        return

    df = pd.DataFrame(rows)
    df = df.sort_values('alignment_score', ascending=False).reset_index(drop=True)
    df.index = df.index + 1
    df.index.name = 'rank'

    cols = [
        'file', 'num_fragments', 'fragment_lengths', 'gap_lengths',
        'fragment_ranges', 'ref_range', 'identity_percent',
        'query_coverage_percent', 'alignment_score'
    ]
    df = df[cols]

    df.to_csv(output, sep="\t", index=True)
    print(f"\nSaved to {output}")
    print(df.head(15).to_string())


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: python script.py <reference.fasta> <pdb_directory/> <output.tsv>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
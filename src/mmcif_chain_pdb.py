#!/usr/bin/env python3
"""
Extract a single chain from an mmCIF file and save it as PDB-format
atomic coordinates (ATOM/HETATM records) to --output.

Usage:
    python mmcif_chain_pdb.py --input path/to/structure.cif --chain_id A --output path/to/out.pdb

Requires:
    pip install biopython
"""

#### imports
import argparse
import io
import os
import sys
from Bio.PDB import MMCIFParser, PDBIO, Select


#### hardcoded values
QUIET = True


#### classes
class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.get_id() == self.chain_id


#### functions
def get_chain_pdb_text(cif_path: str, chain_id: str) -> str:
    parser = MMCIFParser(QUIET=QUIET)
    structure = parser.get_structure("structure", cif_path)

    found = any(chain_id in model for model in structure)
    if not found:
        raise ValueError(f"Chain '{chain_id}' not found in {cif_path}")

    io_writer = PDBIO()
    io_writer.set_structure(structure)

    buffer = io.StringIO()
    io_writer.save(buffer, select=ChainSelect(chain_id))
    return buffer.getvalue()


#### CLI
def main():
    parser = argparse.ArgumentParser(
        description="Extract a chain from mmCIF and save it in PDB format."
    )
    parser.add_argument("--input", required=True, help="Path to input mmCIF file")
    parser.add_argument("--chain_id", required=True, help="Chain ID to extract")
    parser.add_argument("--output", required=True, help="Full output file path")
    args = parser.parse_args()

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    pdb_text = get_chain_pdb_text(args.input, args.chain_id)

    with open(args.output, "w") as fh:
        fh.write(pdb_text)

    print(args.output)


if __name__ == "__main__":
    main()

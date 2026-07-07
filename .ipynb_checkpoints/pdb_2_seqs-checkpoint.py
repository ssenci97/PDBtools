import os
import pandas as pd
from pathlib import Path
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

def get_pdb_sequence(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    
    sequence_parts = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " ":
                    try:
                        sequence_parts.append(seq1(residue.get_resname()))
                    except Exception:
                        continue
                        
    return "".join(sequence_parts)

def main():
    base_dir = Path("data/pdbs_query/chains")
    pdb_files = list(base_dir.glob("*.pdb"))
    
    data = []
    
    for pdb_file in pdb_files:
        fname = pdb_file.stem                
        try:
            sequence = get_pdb_sequence(str(pdb_file))
            
            parts = fname.split("_")
            complex_name = parts[0] if len(parts) > 0 else fname
            pdb_id = parts[1] if len(parts) > 1 else fname
            
            data.append({
                "fname": fname,
                "complex": complex_name,
                "id": pdb_id,
                "sequence": sequence
            })
        except Exception as e:
            print(f"Error processing {pdb_file}: {e}")

    output_path = "data/seqs_pdbs/inference_input_seqs.csv"
    output_dir = Path(output_path).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False)
    print(f"Successfully saved dataframe to {output_path}")

    for entry in data:
        fasta_name = f"{entry['fname']}.fasta"
        fasta_path = output_dir / fasta_name
        
        with open(fasta_path, "w") as f:
            f.write(f">{fasta_name.replace('.fasta', '')}\n")
            f.write(f"{entry['sequence']}\n")

if __name__ == "__main__":
    main()
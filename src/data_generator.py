import random

def generate_synthetic_fasta(filename, num_seqs=10, seq_len=150, motif="GATTACA"):
    """
    Generates synthetic DNA sequences and plants a known motif in each.
    """
    bases = ['A', 'C', 'G', 'T']
    motif_len = len(motif)
    true_positions = []

    with open(filename, 'w') as f:
        for i in range(num_seqs):
            # Generate random background sequence
            seq_list = random.choices(bases, k=seq_len)
            
            # Plant the motif at a random position
            insert_pos = random.randint(0, seq_len - motif_len)
            seq_list[insert_pos:insert_pos + motif_len] = list(motif)
            
            # Convert back to string
            final_seq = "".join(seq_list)
            true_positions.append(insert_pos)
            
            # Write to FASTA format
            f.write(f">seq_{i} planted_at_{insert_pos}\n")
            f.write(f"{final_seq}\n")

    print(f"Generated {num_seqs} sequences with motif '{motif}' saved to {filename}")
    print(f"Ground Truth Positions: {true_positions}")

if __name__ == "__main__":
    # Generate the prototype dataset
    generate_synthetic_fasta("../data/synthetic_promoters.fasta")

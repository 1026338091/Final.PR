import random
import numpy as np

def read_fasta(filename):
    """read FASTA file，extract 10 seq"""
    seqs = []
    with open(filename, 'r') as f:
        seq = ""
        for line in f:
            if line.startswith(">"):
                if seq: seqs.append(seq)
                seq = ""
            else:
                seq += line.strip()
        if seq: seqs.append(seq)
    return seqs

def build_pwm(motifs, k, pseudocount=1):
    """
   setp 3 Construct the position weight matrix PWM
    introduce pseudocount  It's like adding a background noise to the laboratory equipment to prevent the program from crashing due to a zero probability calculation
    """
    counts = {'A': np.ones(k) * pseudocount,
              'C': np.ones(k) * pseudocount,
              'G': np.ones(k) * pseudocount,
              'T': np.ones(k) * pseudocount}
    
    num_motifs = len(motifs)
    for motif in motifs:
        for i, nucleotide in enumerate(motif):
            counts[nucleotide][i] += 1
            
    # Transform into a probability distribution
    total_counts = num_motifs + (4 * pseudocount)
    pwm = {base: counts[base] / total_counts for base in ['A', 'C', 'G', 'T']}
    return pwm

def score_kmer(kmer, pwm):
    """step 4：Score with the standard (calculate the probability of k-mer under PWM)"""
    prob = 1.0
    for i, nucleotide in enumerate(kmer):
        prob *= pwm[nucleotide][i]
    return prob

def roulette_wheel_selection(probabilities):
    """step 5:Roulette selection (Randomly drawn based on probability, the plots with higher scores are selected for larger areas)"""
    total = sum(probabilities)
    if total == 0:
        return random.randint(0, len(probabilities) - 1) # in case erro
        
    normalized_probs = [p / total for p in probabilities]
    return np.random.choice(len(probabilities), p=normalized_probs)

def gibbs_sampler(seqs, k, iterations=150):
    """grand cycle"""
    t = len(seqs)
    n = len(seqs[0])
    
    # Blind selection with closed eyes: Randomly pick a k-mer from each sequence as the initial candidate
    motifs = []
    for seq in seqs:
        start = random.randint(0, n - k)
        motifs.append(seq[start:start+k])
        
    for i in range(iterations):
        # One retention method: Randomly cover a field (eliminate one sequence)
        remove_idx = random.randint(0, t - 1)
        left_out_seq = seqs[remove_idx]
        remaining_motifs = motifs[:remove_idx] + motifs[remove_idx+1:]
        
        # Analyze the remaining fields to obtain the current characteristics (PWM)
        pwm = build_pwm(remaining_motifs, k)
        
        # Go to the covered field and score every inch of land
        probs = []
        for j in range(n - k + 1):
            kmer = left_out_seq[j:j+k]
            prob = score_kmer(kmer, pwm)
            probs.append(prob)
            
        # Roulette draws lots to select new candidate plots and replace the original ones
        new_start = roulette_wheel_selection(probs)
        motifs[remove_idx] = left_out_seq[new_start:new_start+k]
        
    return motifs

if __name__ == "__main__":
    print(" act Gibbs Sampler prototype...")
    # the length of GATTACA we find is 7
    k_len = 7 
    try:
        # Make sure the path points to the test data you generated earlier
        sequences = read_fasta("../data/synthetic_promoters.fasta")
        
        # Iterate 150 times to allow the algorithm to self-correct
        final_motifs = gibbs_sampler(sequences, k=k_len, iterations=150)
        
        print("\n the final filed：")
        for idx, m in enumerate(final_motifs):
            print(f"filed {idx}: {m}")
            
    except FileNotFoundError:
        print(" can't find data,runing data_generator.py first！")

"""
gibbs.py
The main entry point for the MCMC Gibbs Sampler algorithm.
Integrates Bi-directional scanning and Log-Odds Background penalties.
"""
import random
import argparse
import yaml
import numpy as np
from sequence_utils import extract_kmers_both_strands
from background import train_background_model, get_background_prob

def read_fasta(filename):
    """Parses a multi-FASTA file."""
    seqs = []
    with open(filename, 'r') as f:
        seq = ""
        for line in f:
            if line.startswith(">"):
                if seq: seqs.append(seq.upper())
                seq = ""
            else:
                seq += line.strip().upper()
        if seq: seqs.append(seq.upper())
    return seqs

def build_pwm(motifs, k, pseudocount=0.1):
    """Constructs a Position Weight Matrix (PWM) from a list of motifs."""
    pwm = {b: np.ones(k) * pseudocount for b in ['A', 'C', 'G', 'T']}
    total = len(motifs) + (4 * pseudocount)
    
    for motif in motifs:
        for i, base in enumerate(motif):
            if base in pwm:
                pwm[base][i] += 1
                
    for base in pwm:
        pwm[base] /= total
    return pwm

def score_pwm_prob(kmer, pwm):
    """Calculates P(kmer | PWM)."""
    prob = 1.0
    for i, base in enumerate(kmer):
        if base in pwm:
            prob *= pwm[base][i]
        else:
            prob *= 0.25 
    return prob

def run_single_sampler(seqs, k, bg_model, bg_order, iterations, burn_in):
    """Executes a single chain of the Gibbs Sampler."""
    t = len(seqs)
    
    motifs = []
    for seq in seqs:
        start = random.randint(0, len(seq) - k)
        motifs.append(seq[start:start+k])
        
    best_motifs = motifs.copy()
    best_score = -float('inf')
        
    for i in range(iterations):
        remove_idx = random.randint(0, t - 1)
        left_out_seq = seqs[remove_idx]
        remaining_motifs = motifs[:remove_idx] + motifs[remove_idx+1:]
        
        pwm = build_pwm(remaining_motifs, k)
        
        candidates = extract_kmers_both_strands(left_out_seq, k)
        weights = []
        
        for kmer, start_idx, strand in candidates:
            p_pwm = score_pwm_prob(kmer, pwm)
            p_bg = get_background_prob(kmer, bg_model, bg_order)
            weight = p_pwm / p_bg if p_bg > 0 else 0
            weights.append(weight)
            
        total_weight = sum(weights)
        if total_weight == 0:
            selected_idx = random.randint(0, len(candidates) - 1)
        else:
            normalized = [w / total_weight for w in weights]
            selected_idx = np.random.choice(len(candidates), p=normalized)
            
        new_motif = candidates[selected_idx][0]
        motifs[remove_idx] = new_motif
        
        if i >= burn_in:
            current_log_likelihood = sum(np.log2(score_pwm_prob(m, pwm)) for m in motifs)
            if current_log_likelihood > best_score:
                best_score = current_log_likelihood
                best_motifs = motifs.copy()
                
    return best_motifs, best_score

def main():
    parser = argparse.ArgumentParser(description="MCMC Gibbs Sampler for Motif Discovery")
    parser.add_argument("--config", required=True, help="Path to config YAML file")
    parser.add_argument("--data", required=True, help="Path to input FASTA file")
    args = parser.parse_args()
    
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)['algorithm']
        
    k = config['k_width']
    order = config['background_order']
    restarts = config['restarts']
    
    print(f"Initializing Sampler on {args.data}")
    print(f"Parameters: k={k}, order={order}, restarts={restarts}")
    
    sequences = read_fasta(args.data)
    if not sequences:
        print("Error: No sequences found in FASTA file.")
        return
        
    print("Training Higher-Order Markov Background Model...")
    bg_model = train_background_model(sequences, order=order)
    
    global_best_score = -float('inf')
    global_best_motifs = []
    
    print("Starting Random Restarts...")
    for r in range(restarts):
        motifs, score = run_single_sampler(
            seqs=sequences, 
            k=k, 
            bg_model=bg_model, 
            bg_order=order, 
            iterations=config['iterations'], 
            burn_in=config['burn_in']
        )
        print(f"   Restart {r+1}/{restarts} | Log-Likelihood Score: {score:.2f}")
        
        if score > global_best_score:
            global_best_score = score
            global_best_motifs = motifs
            
    print("\nProcess Complete. Best Motif Instances Discovered:")
    for i, m in enumerate(global_best_motifs):
        print(f"Seq {i}: {m}")
        
if __name__ == "__main__":
    main()

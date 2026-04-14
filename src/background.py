"""
background.py
Implements the Higher-Order Markov Background Model to dynamically calculate 
compositional bias directly from input data.
"""
from collections import defaultdict

def train_background_model(sequences, order=1, pseudocount=1.0):
    """
    Trains a Markov model of the specified order based on input sequences.
    """
    counts = defaultdict(lambda: defaultdict(float))
    context_totals = defaultdict(float)
    
    for seq in sequences:
        for i in range(len(seq) - order):
            context = seq[i:i+order] if order > 0 else ""
            next_base = seq[i+order]
            counts[context][next_base] += 1
            context_totals[context] += 1
            
    probabilities = defaultdict(dict)
    bases = ['A', 'C', 'G', 'T']
    
    for context in counts.keys() if order > 0 else [""]:
        for base in bases:
            numerator = counts[context].get(base, 0.0) + pseudocount
            denominator = context_totals[context] + (4.0 * pseudocount)
            probabilities[context][base] = numerator / denominator
            
    if order > 0:、
      
        probabilities["_base_"] = base_probs[""]
        
    return probabilities

def get_background_prob(kmer, model, order):
    """
    Calculates the probability of observing a k-mer given the background model.
    """
    prob = 1.0
    if order == 0:
        for base in kmer:
            prob *= model[""].get(base, 0.25)
    elif order == 1:
        prob *= model["_base_"].get(kmer[0], 0.25)
        for i in range(len(kmer) - 1):
            context = kmer[i]
            next_base = kmer[i+1]
            prob *= model.get(context, {}).get(next_base, 0.25)
            
    return prob

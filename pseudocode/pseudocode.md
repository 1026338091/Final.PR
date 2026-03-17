
# Part 2: Conceptual Progress Report (Pseudocode & Design) Hongyuan Deng

1. Conceptual Overview
This project implements a **De Novo Motif Discovery tool** using the Gibbs Sampling algorithm (a Markov Chain Monte Carlo approach). To address the common pitfall of stochastic algorithms getting trapped in local optima, the design incorporates a **Random Restart** wrapper. Furthermore, to handle the compositional bias of real eukaryotic promoters (e.g., GC-rich regions), the core scoring mechanism evaluates candidate k-mers against a **Background Model** (log-odds scoring) rather than assuming a uniform nucleotide distribution.

---

 2. Detailed Pseudocode

 Module 2.1: Main Wrapper (Handling Local Optima via Restarts)
This outer loop executes the core Gibbs Sampler multiple times with different random initializations to ensure robust global convergence.

    FUNCTION GibbsMotifFinder_WithRestarts(DnaSequences, k, num_restarts, max_iterations):
        INPUT: 
            DnaSequences: List of N promoter strings
            k: Integer, length of the motif
            num_restarts: Integer, number of independent runs
            max_iterations: Integer, max steps per run to reach convergence
        OUTPUT: 
            best_overall_motifs: List of the best k-mers found
            best_overall_pwm: The corresponding Position Weight Matrix

        best_overall_score = -INFINITY
        best_overall_motifs = NULL
        best_overall_pwm = NULL

        FOR run FROM 1 TO num_restarts:
            // Execute single Gibbs sampling run
            current_motifs, current_pwm, current_score = RunSingleGibbs(DnaSequences, k, max_iterations)

            // Update global best if the current run converged to a higher-scoring state
            IF current_score > best_overall_score:
                best_overall_score = current_score
                best_overall_motifs = current_motifs
                best_overall_pwm = current_pwm

        RETURN best_overall_motifs, best_overall_pwm


Module 2.2: Core Gibbs Sampler (Leave-One-Out Optimization)
The stochastic heart of the algorithm. It iteratively refines the motif by leaving one sequence out, building a model from the rest, and probabilistically sampling a new position in the left-out sequence.

    FUNCTION RunSingleGibbs(DnaSequences, k, max_iterations):
        N = length(DnaSequences)
        
        // 1. Initialization
        motifs = RandomlySelectInitialKmers(DnaSequences, k)
        
        // 2. Background Model Initialization (Key Innovation)
        // Calculates single-nucleotide frequencies (0-order) or dinucleotide frequencies across all input data
        bg_model = TrainBackgroundModel(DnaSequences) 

        best_motifs = motifs
        best_score = CalculateInformationContent(motifs, bg_model)

        FOR iter FROM 1 TO max_iterations:
            // Select a sequence to hold out randomly
            j = RandomInteger(1, N)
            left_out_seq = DnaSequences[j]
            
            // Build PWM from the (N-1) remaining sequences
            // Note: Pseudocounts are added to avoid zero probabilities
            pwm = BuildPWM_WithPseudocounts(motifs EXCLUDING motifs[j])

            // Score all possible starting positions in the left-out sequence
            probabilities = EMPTY LIST
            FOR EACH valid starting position 'pos' IN left_out_seq:
                kmer = ExtractKmer(left_out_seq, pos, k)
                
                // Calculate Log-Odds score: log( P(kmer|PWM) / P(kmer|Background) )
                prob_motif = CalculateProbability(kmer, pwm)
                prob_bg = CalculateProbability(kmer, bg_model)
                
                weight = prob_motif / prob_bg
                APPEND weight TO probabilities

            // Normalize weights into a true probability distribution (sum = 1.0)
            normalized_probs = Normalize(probabilities)

            // Probabilistic Sampling (Roulette Wheel Selection)
            // We do NOT pick the maximum value; we sample based on the probability distribution
            new_pos = RandomChoiceBasedOnProbabilities(normalized_probs)
            motifs[j] = ExtractKmer(left_out_seq, new_pos, k)

            // Evaluate convergence
            current_score = CalculateInformationContent(motifs, bg_model)
            IF current_score > best_score:
                best_score = current_score
                best_motifs = motifs
            ELSE IF score_has_not_improved_for_X_iterations:
                BREAK // Early stopping

        RETURN best_motifs, BuildPWM(best_motifs), best_score

---

 3. Complexity Analysis

Time Complexity: In a single iteration, building the PWM takes O(N * k) operations. Scanning the left-out sequence of length L to score all (L - k + 1) possible k-mers takes O(L * k) operations. For a single Gibbs run with I iterations, the complexity is O(I * L * k). Factoring in R random restarts, the overall time complexity is O(R * I * L * k). 
  *Tradeoff discussion:* Since k (motif length, typically 6-15) and R (restarts, typically 20-50) are small constants, the algorithm scales linearly with the sequence length L and the number of iterations I. This makes it computationally highly efficient and feasible to run locally without needing an HPC cluster for standard promoter sets.

pace Complexity: The algorithm requires memory to store the sequences O(N * L), the current set of motifs O(N * k), the PWM O(4 * k), and the background model. The overall space complexity is tightly bounded by the input size: O(N * L).
  *Tradeoff discussion:* Memory is not a bottleneck for this algorithm. However, moving to higher-order Markov models (e.g., 3rd-order) would exponentially increase the space needed for the background probability table, which is a manageable but necessary tradeoff for biological accuracy.

---

4. Testing and Validation Plan

Before applying the algorithm to messy biological data, I will build a progressive testing pipeline to ensure the core logic is sound.

 Phase 1: Unit Testing (Sanity Checks)
 PWM Construction Check: Pass a known set of identical k-mers (e.g., five `ATGCGT` sequences) to the `BuildPWM` function and verify that the resulting matrix outputs probabilities of 1.0 for the consensus bases and accounts for pseudocounts correctly.
Sampling Check: Force the probability distribution array to `[0.0, 1.0, 0.0, 0.0]`. The `RandomChoice` function must deterministically return index 1 exactly 100% of the time.

Phase 2: Integration Testing (Synthetic Ground Truth)
Action: I will write a script to generate 10 random DNA sequences (length 150bp) based on a uniform distribution. I will then randomly "plant" a specific motif (e.g., `GATTACA`) into each sequence, allowing for 1-2 random point mutations to simulate biological variation.
Success Metric: The algorithm must recover the starting positions of the planted motifs with >= 90% accuracy. 

Phase 3: Edge Case Testing
Phase Shift Tolerance: I will write an evaluation function that counts a predicted motif as "correct" if its start index is within +/- 2 bp of the true planted index, acknowledging that Gibbs sampling often converges on slightly shifted overlapping sequences.
Adversarial Background: I will generate synthetic sequences with a heavy GC bias (e.g., 75% G/C) and plant an AT-rich motif. This will explicitly test if my Background Model correctly penalizes GC-rich noise and successfully isolates the true motif.

###AI Usage：Use Genmini and Claud for Format Building

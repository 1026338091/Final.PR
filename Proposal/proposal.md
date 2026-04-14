[proposal.md](https://github.com/user-attachments/files/25377887/proposal.md)
# Project Proposal: From beginning Motif Discovery in Eukaryotic Promoters using Gibbs Sampling with Higher-Order Background Models

## Research Question
**Question:** How can we accurately identify transcription factor binding sites (motifs) in the noisy, non-uniform upstream promoter regions of co-expressed genes, and how does accounting for multi-order genomic background bias improve discovery rates?

**Specific Scope & Biological Rationale:**
To make this research grounded and testable, I will focus on two specific co-expressed gene clusters in *Saccharomyces cerevisiae* (Yeast), targeting approximately **10-15 genes per cluster**
**The Gcn4 Cluster (Amino Acid Starvation):** Gcn4 is a master transcriptional regulator with a well-characterized, relatively neutral motif (`TGACTC`). This 10-15 gene cluster will serve as my **"Gold Standard" positive control** to verify the baseline functionality of the sampler.
**The Mbp1 Cluster (Cell Cycle G1/S Transition):** Mbp1 binds to a highly **GC-rich motif** (`ACGCGT`). Analyzing the promoters of this specific cluster serves as a **"Stress Test"** for the algorithm. It directly addresses the biological rationale of this project: demonstrating whether the proposed Markov background model can successfully prevent the algorithm from being confused by general GC-repeat noise common in eukaryotic genomes.

**Significance & Innovation:** Identifying regulatory elements is key to understanding gene networks, but experimental determination is costly. While standard Gibbs Sampling works well on synthetic data, it often fails in real sequences due to compositional bias. By explicitly testing on the Mbp1 cluster, this project will innovatively demonstrate how integrating a **Higher-Order Markov Background Model** (penalizing candidate sequences based on background transition probabilities) significantly reduces false-positive discoveries compared to naive uniform models.


## Algorithm and Algorithm Class
* **Algorithm Class:** Randomized / Combinatorial Algorithms (Markov Chain Monte Carlo - MCMC)
* **Specific Algorithm:** Gibbs Sampler for Motif Finding (with Background Model correction)
* **Justification & Integration Mechanism:** The motif finding problem is NP-hard. While Gibbs Sampling is a preferred heuristic for its balance of speed and accuracy, its naive form is highly sensitive to genomic noise. To address this, I will integrate a progressively higher-order Markov Chain background model directly into the sampling mechanism:
  1. **Estimation:** The background probability model (e.g., dinucleotide transition frequencies for a 1st-order model) will be estimated directly from the **input promoter sequences** prior to the sampling iterations. This ensures the model dynamically captures the specific, local compositional bias of the analyzed gene cluster.
  2. **Scoring Modification:** During the leave-one-out step, a standard Gibbs sampler scores a candidate DNA sequence (a "k-mer") strictly by calculating how well it matches the current motif model (the Position Weight Matrix, or PWM). My implementation modifies this by dividing that basic "Motif Match Score" by a "Background Probability Score" calculated from the Markov chain. This creates a Log-Odds probability ratio.
  3. **Sampling Integration:** This **Log-Odds ratio** will be used *actively during the sampling phase* to define the probability weights for the roulette wheel selection (not just as a post-hoc evaluation metric). By directly penalizing candidate k-mers that mimic the background noise during every iteration, the sampler is mathematically steered away from ubiquitous repeats and forced to converge on functionally unique regulatory elements.

## Data Plan
### 1. Data Source & Type
**Real Data:** I will extract **full upstream promoter sequences** (e.g., -1000bp to +200bp relative to the Transcription Start Site) for specific co-expressed gene clusters.
**Biological & Methodological Justification:**
    **Why Yeast (*Saccharomyces cerevisiae*)?** Yeast is a premier model eukaryote. It possesses a compact, extensively annotated genome with a highly mapped transcriptional regulatory network. It provides the perfect balance for this project: it is biologically complex enough to contain real genomic noise, yet computationally tractable for a local machine. Most importantly, it has abundant, experimentally validated PWMs available in databases like JASPAR for rigorous ground-truth comparison.
    **Why Promoters instead of ChIP-seq?** Pre-aligned ChIP-seq data artificially simplifies the motif discovery problem by narrowing the search space to pre-enriched fragments centered tightly on the binding event. By using full upstream promoter sequences, I am simulating the true *de novo* biological search problem. The larger search space introduces realistic compositional noise, which is absolutely necessary to test the core innovation of this project: the Higher-Order Background Model.
    **Why Gcn4 and Mbp1 specifically?** As outlined in the Research Question, Gcn4 provides a well-characterized, relatively neutral motif serving as a positive control baseline. Mbp1 binds a highly GC-rich motif, serving as the necessary "stress test" to explicitly demonstrate the background model's ability to filter out GC-repeat noise.
**Format:** Multi-FASTA format, retrieved via YeastMine or SCPD.

### 2. Prototype Data (for initial testing)
To validate the logic before handling biological noise, I will generate **Synthetic Data** with parameters that strictly mirror the expected real biological datasets to avoid scale-up failures:
* **Content:** 10-20 random sequences with lengths matching actual yeast promoters (**~800-1000bp**) with "planted" motifs.
* **Complexity Levels:**
    1.  *Basic:* Uniform background, strong motif.
    2.  *Adversarial:* Highly biased backgrounds designed to trick naive algorithms. This will include not only global GC-bias (e.g., 70% GC content) but also localized low-complexity regions such as **poly-A/T tracts**. (Note: While extreme complexities like nucleosome positioning will not be synthetically simulated, the adversarial data will be rigorous enough to validate the Markov background model before transitioning to real *in vivo* noise).

## Success Criteria
1.  **Deliverable:** A Python command-line tool (CLI) that accepts a FASTA file, motif width $k$, and background model order $m$ (0-order vs. 1st-order, etc.).
2.  **Robustness & Convergence:** To mathematically ensure the MCMC algorithm escapes local optima, the tool must execute and compare a minimum of **50 independent random restarts** per dataset, reporting the consensus motif with the highest global log-odds score.
3.  **Validation:**
    * **Synthetic:** Recover planted motifs with >95% accuracy in both basic and 1000bp adversarial datasets.
    * **Biological:** The tool must identify the known consensus sequences for Gcn4 and Mbp1. The resulting Position Weight Matrices (PWMs) will be aligned and compared against the JASPAR database.
4.  **Visualization:** A dynamic plot showing the increase in Information Content (IC) per iteration, demonstrating algorithm convergence.

## Pitfall Scan

### 1. Algorithmic Issue: Local Optima & MCMC Convergence
* **Concern:** As a stochastic MCMC method, a naive Gibbs Sampler is highly sensitive to initialization. It risks getting trapped in local optima, collapsing into a single strand orientation, or functioning as a greedy hill-climber if the sampling step is implemented deterministically.
* **Mitigation:** To ensure true stochastic exploration and robust convergence, my implementation employs a three-pronged strategy:
    1.  **True Probabilistic Sampling:** During the leave-one-out step, the algorithm will *not* select the maximum-scoring k-mer (i.e., avoiding `argmax`). Instead, it calculates the posterior probability distribution for all valid positions and utilizes **Roulette Wheel Selection** to sample the next motif. This prevents the algorithm from becoming a greedy, easily trapped hill-climber.
    2.  **Bi-directional Strand Evaluation:** Transcription factors bind to double-stranded DNA. At each iteration, the sampler will evaluate candidate k-mers on both the **forward and reverse-complement strands**. The log-odds scores from both strands will be integrated into the single probability distribution to ensure the sampler does not collapse into a single orientation.
    3.  **Burn-in Period & Random Restarts:** To address the instability of early MCMC iterations, each sampling run will include a **Burn-in period** (e.g., ignoring the first 50 iterations) before recording motif statistics or assessing convergence. Furthermore, the entire process will be wrapped in a minimum of **50 independent Random Restarts** to escape local optima and guarantee global convergence.

### 2. Data Issue: Universal Compositional Bias & Circularity
* **Concern:** Eukaryotic promoters contain ubiquitous structural biases beyond just GC-richness, including poly-A/T tracts and low-complexity regions. A naive sampler will falsely lock onto these. Furthermore, estimating the background model directly from the input sequences raises a valid concern of "circularity" (where the true motif slightly skews the background frequencies).
* **Mitigation:**
    * **Addressing Circularity:** In a typical 1000bp yeast promoter, a 10bp motif represents $<1\%$ of the sequence space. Thus, estimating the Markov transition matrix from the entire input cluster is statistically robust and dynamically captures the *local* structural bias. To rigorously control for any residual circularity, I will cross-validate the results by running the sampler using a pre-calculated, genome-wide yeast promoter background model.
    * **Addressing Bias (The Penalty Mechanism):** The Markov model is not specific to GC-bias; it captures *all* transition probabilities (e.g., the high probability of 'A' following 'A' in poly-A tracts). This is actively incorporated during the sampling phase, not post-hoc. The roulette wheel selection utilizes a **Log-Odds Score**: 
      $W = \frac{P(\text{k-mer} \mid \text{PWM})}{P(\text{k-mer} \mid \text{Background})}$
    * By placing the background probability in the denominator, any sequence that scores high simply because it matches the ubiquitous background structure (whether GC-rich or AT-rich) will have its sampling weight heavily penalized, mathematically forcing the algorithm to seek true, functional motifs.

### 3. Evaluation Issue: Phase Shifts & Strand Orientation
* **Concern:** The algorithm might converge on a motif that is shifted by 1-2 bases (finding `ACTCA` instead of `GACTCA`), or it might successfully identify the valid reverse-complement of the true motif. Standard naive metrics checking for exact forward-strand indices will falsely mark these biologically valid recoveries as failures.
* **Mitigation:**
    * **Phase Shift Tolerance:** I will implement a sliding-window overlap check in the evaluation script. A predicted motif will be counted as a "success" if the overlapping aligned region constitutes **>= 80% of the total motif width (k)** (a shift of up to 2 bases for a 10bp motif).
    * **Strand-Agnostic Evaluation:** Aligning directly with the bi-directional sampling approach outlined above, the evaluation metric will explicitly score both the predicted motif and its **reverse-complement** against the ground truth (whether a planted synthetic motif or a JASPAR PWM). A high-scoring overlap in either orientation will be correctly recorded as a true positive, ensuring the algorithm's performance is not artificially penalized.

## Planned Repository Structure (Initial Sketch)
Based on constructive peer review, the repository architecture has been expanded to strictly separate configuration, core logic, visualization, and logging, ensuring reproducibility and easier debugging of the MCMC process.

```text
/
├── README.md               # Project overview, installation, and Quick Start
├── PROPOSAL.md             # This document
├── environment.yml         # Conda environment dependencies
├── configs/                # NEW: Configuration files (YAML) for hyperparameters
│   └── default_config.yml  # Stores motif width (k), max iterations, restarts, etc.
├── src/                    # Source code package
│   ├── gibbs.py            # Main Gibbs Sampling logic with Random Restarts
│   ├── background.py       # Markov Chain probability matrix calculations
│   ├── scoring.py          # Log-odds probability ratio calculations
│   ├── sequence_utils.py   # NEW: Bi-directional scanning & reverse-complement logic
│   ├── synthetic.py        # NEW: Programmatic generation of prototype/adversarial data
│   ├── visualize.py        # NEW: Modular plotting for IC convergence and PWM logos
│   └── io_utils.py         # FASTA parsing and CSV logging
├── data/
│   ├── raw/                # Yeast promoter FASTA files
│   └── benchmarks/         # Synthetic datasets generated by synthetic.py
├── models/                 # NEW: Saved Markov background transition matrices
├── outputs/                # NEW: Tracked MCMC logs, intermediate PWMs, and IC history
├── tests/                  # Unit tests
│   ├── test_sampling.py
│   └── test_background.py
└── notebooks/              # High-level analysis importing cleanly from src/



Generative AI Disclosure
Tool Used: Gemini

Drafting: Used to articulate the Pitfall Scan, specifically regarding genomic background bias and phase shifts.

Differentiation Strategy: Used to brainstorm ways to distinguish this project from standard assignments by introducing Higher-Order Markov Models.

Formatting: Use it to adjust the structure of my proposal to make it more professional

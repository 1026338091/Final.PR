[proposal.md](https://github.com/user-attachments/files/25377887/proposal.md)
# Project Proposal: From beginning Motif Discovery in Eukaryotic Promoters using Gibbs Sampling with Higher-Order Background Models

#Hongyuan Deng

## Research Question
**Question:** How can we accurately identify transcription factor binding sites (motifs) in the noisy, non-uniform upstream promoter regions of co-expressed genes, and how does accounting for genomic background bias improve discovery rates?

**Significance:** Identifying regulatory elements is key to understanding gene networks. While standard Gibbs Sampling works well on synthetic data, it often fails in real genomic sequences due to compositional bias (e.g., GC-rich regions, CpG islands).
**Innovation:** Unlike basic implementations that assume a uniform background distribution, this project will implement a **Gibbs Sampler with a Higher-Order Markov Background Model**. I will systematically compare how a 0-order vs. a 3rd-order (dinucleotide/trinucleotide)background model affects the ability to filter out false positives in complex eukaryotic promoters.

## Algorithm and Algorithm Class
* **Algorithm Class:** Randomized / Combinatorial Algorithms (Markov Chain Monte Carlo - MCMC)
* **Specific Algorithm:** Gibbs Sampler for Motif Finding (with Background Model correction)
* **Justification:** The motif finding problem is NP-hard. Gibbs Sampling is a preferred heuristic for its balance of speed and accuracy. However, its standard form is sensitive to genomic noise. By integrating a Markov chain to model the "background" probability of a sequence, the algorithm can distinguish true motifs from random repetitive elements more effectively.

## Data Plan
### 1. Data Source & Type
* **Real Data:** I will not use pre-aligned ChIP-seq reads. Instead, I will extract **full upstream promoter sequences** for specific gene clusters in *Saccharomyces cerevisiae* (Yeast).
    * **Source:** The SCPD(Promoter Database of Saccharomyces cerevisiae) or extracted via UCSC Genome Browser.
    * **Target:** Gene sets regulated by **Gcn4** (well-characterized) and **Mbp1** (cell cycle regulator).
    * **Format:** FASTA format.

### 2. Prototype Data
To validate the logic before handling biological noise, I will generate **Synthetic Data**:
* **Content:** 10-20 random sequences (length 200bp) with "planted" motifs.
* **Complexity Levels:**
    1.  Basic: Uniform background, strong motif.
    2.  Advanced: Biased background (e.g., 70% GC content) to test if the algorithm gets confused without the Markov model.

## Success Criteria
1.  **Deliverable:** A Python command-line tool (CLI) that accepts a FASTA file, motif width $k$, and background model order $m$.
2.  **Validation:**
    * **Synthetic:** Recover planted motifs with >95% accuracy in the "Basic" scenario.
    * **Biological:** The tool must identify the known consensus sequence for Gcn4 (`TGACTC`) from real promoter data, and the Resulting PWM (Position Weight Matrix) will be compared against the JASPAR database using Euclidean distance or Pearson correlation.
3.  **Visualization:** A dynamic plot showing the increase in **Information Content (IC)** per iteration, demonstrating the algorithm's convergence behavior.

## Pitfall Scan

### 1. Algorithmic Issue: Local Optima (The "Stuck" Problem)
* **Concern:** Gibbs Sampling is a stochastic hill-climbing method. It can easily get trapped in a local optimum (a shifted or weak motif) and never find the global best motif.
* **Mitigation:** I will implement **Random Restarts** (running the algorithm multiple times from different random seeds) and report the motif with the highest final Information Content score.

### 2. Data Issue: Genomic Background Bias (The "GC" Problem)
* **Concern:** Real promoters are not random strings of letters. A simple Gibbs Sampler might "discover" a region simply because it is GC-rich, not because it is a motif. This is a major limitation of basic implementations.
* **Mitigation:** This is the core innovation of my project. I will calculate a **Background Markov Model** from the input sequences. When scoring a potential motif, I will calculate a log-odds score: $log(P(sequence|motif) / P(sequence|background))$, effectively penalizing sequences that just look like the genomic background.

### 3. Evaluation Issue: Phase Shifts (The "Drift" Problem)
* **Concern:** The algorithm might converge on a motif that is shifted by 1-2 bases (e.g., finding `ACTCA` instead of `GACTCA`). Standard accuracy metrics checking for exact indices will mark this as a failure.
* **Mitigation:** I will implement a "Phase Shift" check in my evaluation code. If the discovered motif strongly overlaps the true motif (e.g., overlap > 80%), it will be counted as a success, and I will implement a step to realign the PWM columns if a shift is detected.

## Planned Repository Structure (Initial Sketch)
```text
/
├── README.md               # Project overview, installation, and Quick Start
├── PROPOSAL.md             # This document
├── environment.yml         # Conda environment dependencies
├── src/                    # Source code package
│   ├── __init__.py
│   ├── gibbs.py            # Main Gibbs Sampling logic
│   ├── background.py       # Markov Chain background modeling (New)
│   ├── scoring.py          # PWM and Information Content calculations
│   └── io_utils.py         # FASTA parsing
├── data/
│   ├── raw/                # Yeast promoter FASTA files
│   ├── benchmarks/         # Synthetic datasets with known answers
│   └── jaspar_motifs/      # Ground truth PWMs for validation
├── tests/                  # Unit tests
│   ├── test_sampling.py
│   └── test_background.py
└── notebooks/              # Analysis of Background Model impact
    └── comparison.ipynb




Generative AI Disclosure
Tool Used: Gemini
Usage:

Drafting: Used to articulate the Pitfall Scan, specifically regarding genomic background bias and phase shifts.

Differentiation Strategy: Used to brainstorm ways to distinguish this project from standard assignments by introducing Higher-Order Markov Models.

Formatting: Use it to adjust the structure of my proposal to make it more professional

# De Novo Motif Discovery with Higher-Order Markov Background Models

## Project Overview
* **Research Question:** How can we accurately identify transcription factor binding sites (motifs) in the noisy, non-uniform upstream promoter regions of co-expressed genes, and how does accounting for multi-order genomic background bias improve discovery rates?
* **Algorithm Class:** Randomized / Combinatorial Algorithms (Markov Chain Monte Carlo - MCMC).
* **Specific Algorithm:** Gibbs Sampler integrated with a dynamically trained Higher-Order Markov Background Model and Bi-directional strand scanning.
* **Data & Outputs:** Analyzes multi-FASTA files of *Saccharomyces cerevisiae* (Yeast) promoters (e.g., Gcn4 and Mbp1 clusters). The output includes the discovered consensus sequence, the aligned Position Weight Matrix (PWM), and Information Content (IC) progression logs.

## Repository Structure
This repository is organized to separate configuration, core algorithm logic, and data, ensuring reproducibility.

* `/configs/`: YAML configuration files for hyperparameters.
* `/src/`: Core Python modules (Gibbs sampler, background model, strand utilities).
* `/data/`: Input FASTA sequences (synthetic and biological).
* `/tests/`: Automated testing scripts.

## Installation
This tool requires Python 3.9 or higher. To ensure a clean environment, please set up a virtual environment.

**1. Clone the repository and navigate into it:**
git clone https://github.com/1026338091/Final.PR.git
cd Final.PR

**2. Install required dependencies:**
pip install -r requirements.txt

*(Note: The only major external dependencies are `numpy` for probability calculations and `PyYAML` for parsing configuration files).*

## Quick Start
To verify the installation, you can run a minimal end-to-end test using the provided synthetic dataset.

**Command to run:**
python src/gibbs.py --config configs/default.yaml --data data/synthetic.fasta

**Expected Output:**
The script will execute the MCMC sampling process and output the following to the console:
1. Confirmation of the chosen hyperparameters (k-mer length, background order).
2. Training status of the Markov background model.
3. The convergence of Log-Likelihood scores across independent random restarts.
4. The final predicted Motif Consensus and the starting positions of the motif in each sequence.

## Usage and Options
The algorithm's behavior is primarily controlled by the `configs/default.yaml` file to prevent hard-coding. 

**Key Configuration Parameters:**
* `k_width`: Length of the target motif.
* `background_order`: `0` for a standard uniform background, `1` for a dinucleotide Markov model, `2` for trinucleotide.
* `restarts`: Number of independent chains to run to escape local optima.

**Command Line Execution:**
You can override the data file and background order directly from the command line. For example, to analyze real Yeast promoters with a 1st-order model:
python src/gibbs.py --config configs/default.yaml --data data/raw/gcn4_promoters.fasta

## Limitations and Assumptions
* **High-Order Data Sparsity:** A 2nd or 3rd-order Markov model theoretically captures more complex structural noise (like nucleosome positioning). However, on small datasets (e.g., < 10 promoter sequences), these high-order transition matrices will suffer from data sparsity (many zero probabilities), leading to inaccurate log-odds penalties. The tool defaults to a 1st-order model for balance.
* **Fixed Motif Width:** This Gibbs sampler assumes the motif width ($k$) is fixed and known *a priori*. It does not handle motifs with variable-length insertions or deletions (indels).

## Evidence of Correctness
1. **Automated Mathematical Checks:** The `tests/test_background.py` suite explicitly verifies that the generated Markov Transition Probability matrices are valid (i.e., all probability distributions sum exactly to 1.0), preventing probability leaks.
2. **Biological Baseline Validation:** When analyzing the Yeast Gcn4 promoter cluster, the algorithm successfully escapes local optima through random restarts and recovers a consensus sequence that strongly aligns with the established JASPAR database standard (`TGACTC`), demonstrating functional correctness *in vivo*.

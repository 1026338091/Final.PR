# Final Reflection

## What Went Right
One of the biggest successes in this project was successfully pivoting the architecture based on the rigorous peer reviews I received in Part 1 and Part 3. My peers correctly pointed out that a naive Gibbs Sampler would fail against bi-directional DNA properties and hard-coded parameters. By modularizing the code (separating `gibbs.py`, `background.py`, and `sequence_utils.py`) and introducing a `configs/default.yaml` file, the tool evolved from a rigid script into a highly reproducible, professional-grade pipeline. Furthermore, successfully implementing the Log-Odds penalty using the Higher-Order Markov Model proved highly effective in preventing the sampler from locking onto GC-repeat noise during testing.

## What Went Wrong or Was Hard
The mathematical reality of translating probabilistic algorithms into code was far harder than writing the pseudocode. Specifically, I severely underestimated the issue of `log(0)` errors. In early implementations, if a specific k-mer was never seen in the background model, its probability was 0. When calculating the log-likelihood, the program would crash instantly. Implementing proper "pseudocounts" to smooth the data was a major technical hurdle. 
On a personal level, dealing with illness during the critical implementation phase (Part 3) made debugging extremely frustrating. I learned that taking a step back and thinking conceptually (rather than blindly staring at syntax errors) was essential for my mental health and project progress.

## Algorithmic Lessons
This project fundamentally changed my understanding of Markov Chain Monte Carlo (MCMC) methods. In lecture, Gibbs Sampling sounded almost "magical" in its ability to find motifs. In practice, I learned that MCMC algorithms are incredibly "greedy" and sensitive to initialization. Without implementing Bi-directional Strand Scanning and forcing at least 10-20 Random Restarts, the algorithm would constantly get trapped in local optima (e.g., finding a phase-shifted version of the motif). The tradeoff between speed and accuracy is very real: running 50 restarts with a 2nd-order background model yields brilliant results, but the computational time scales up dramatically compared to a naive 0-order run.

## Future Directions
If I were to extend this project, I would implement an automatic sequence-padding feature to allow the discovery of motifs with variable lengths (insertions/deletions), as the current rigid k-mer width is biologically restrictive. Additionally, integrating a visualization library like `matplotlib` or `seaborn` to automatically generate sequence logos (rather than just printing the text consensus) would make the tool much more user-friendly for biologists.

## Generative AI Disclosure
* **Tool Used:** Gemini (Google)
* **Usage Context:** Generative AI was utilized as a technical tutor throughout the project. It was highly instrumental in helping me scaffold the initial `pytest` scripts to verify my probability math, translating my conceptual pseudocode into the final Python syntax (especially the dictionary structures for the Markov transition matrix), and assisting in proofreading and structuring my responses to the peer review feedback in Part 1 to ensure a professional, academic tone.

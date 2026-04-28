# Parallel implementation of global or local sequence alignment (Needlman-Wunsch or Smith-Waterman)

## This project implements a parallelized version of **Global (Needleman-Wunsch)** and **Local (Smith-Waterman)** sequence alignment algorithms using Python's `multiprocessing`.

## Project structure 
- `alignment_script.py`: the main python script
- `example_data/`: folder containing FASTA sequences for the user to try
               -`sequence1.fasta`: human hemoglobin subunit alpha
               -`sequence2.fasta`: human hemoglobin subunit beta

### Features
- *Parallel Matrix Computation*: Exploits multi-core architectures to speed up the scoring matrix filling.
- *Dual Mode*: Support for both global and local alignment.
- *Flexible Scoring*: Customizable match score, mismatch penalty, and gap penalty.
- *FASTA Support*: Reads sequences directly from `.fasta` files.

## Parallelization Strategy

### Theoretical Concept
In dynamic programming for sequence alignment, a cell $(i, j)$ depends on cells $(i-1, j)$, $(i, j-1)$, and $(i-1, j-1)$. This dependency prevents a simple row-by-row or column-by-column parallelization.

However, all cells $(i, j)$ where the sum $i + j$ is constant are independent of each other. These cells form a **minor diagonal**.

### Implementation
The script processes the matrix diagonal by diagonal:
1. It identifies all cells belonging to the current diagonal.
2. It partitions these cells into **chunks** to minimize the overhead of inter-process communication. 
> (technical note): For short sequences, the overhead is negligible. For very long sequences, the chunking strategy is implemented specifically to balance the trade-off between communication overhead and computational speedup.
3. It uses a `ProcessPoolExecutor` to compute the scores of each chunk in parallel across available CPU cores.

**Maximum Parallelism**: The maximum number of cells that can be computed simultaneously is equal to the length of the longest diagonal, which is $\min(N, M)$ (where $N$ and $M$ are sequence lengths). Therefore, using more processes than $\min(N, M)$ does not provide additional speedup.

## Usage
### To test the tool with the provided example data:

1. Clone the repository: 
   ```bash
   git clone [https://github.com/fulviadercole/Parallel-Alignment-Python.git](https://github.com/fulviadercole/Parallel-Alignment-Python)
2. Run the script:
   ```bash
   python alignment_script.py

### Prerequisites
- Python 3.8+
- NumPy (`pip install numpy`)


## The final result will be printed in the console and saved to `alignment_output.txt`. 

## Contact: 
### Fulvia D'Ercole - fulvia.dercole@studenti.unimi.it
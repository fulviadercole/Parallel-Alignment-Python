import numpy as np
from concurrent.futures import ProcessPoolExecutor
import time
import multiprocessing


#input and output functions
def read_fasta(file_path):
    """Reading a sequence from a .fasta file ignoring headers"""
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        sequence = "".join(line.strip() for line in lines if line.strip() and not line.startswith(">"))
        return sequence.upper()
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        return None
    
def save_alignment_to_file(a1, a2, score, mode, filename="alignment_output.txt"):
    with open(filename, "w") as f:
        f.write(f" Alignment result ({mode.upper()}) ---\n")
        f.write(f"Total score: {score}\n")
        f.write("-" * 40 + "\n\n")

        #rapping every 60 characters for readability
        for i in range(0, len(a1), 60):
            f.write(f"SEQ1: {a1[i:i+60]}\n")
            f.write(f"SEQ2: {a2[i:i+60]}\n\n")

    print(f"[*] Result saved successfully in: {filename}")

#Parallel calculation 
def compute_chunk_tasks(args):
    cells_to_compute, seq1, seq2, match, mismatch, gap, matrix, mode = args
    results = []
    for (i, j) in cells_to_compute:
        diag = matrix[i-1, j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
        up = matrix[i-1, j] + gap
        left = matrix[i, j-1] + gap

        val = max(diag, up, left)
        if mode == 'local':
            val = max(0, val)
        results.append((i, j, val))
    return results

def parallel_alignment(seq1, seq2, match=1, mismatch=-1, gap=-1, mode='global'):
    rows, cols = len(seq1) + 1, len(seq2) + 1
    matrix = np.zeros((rows, cols), dtype=np.int32)

    #Edge initialization
    if mode == 'global':
        matrix[:, 0] = np.arange(rows) * gap
        matrix[0, :] = np.arange(cols) * gap

    num_cores = multiprocessing.cpu_count()


    with ProcessPoolExecutor(max_workers=num_cores) as executor:
            for d in range(2, rows + cols):
                #Identifying cells on the current diagonal
                diag_cells = [(i, d-i) for i in range(1, rows) if 1 <= d-i < cols]

                #Chunking to reduce inter-process communication overhead
                chunk_size = max(1, len(diag_cells) // (num_cores * 2))
                chunks = [diag_cells[x:x+chunk_size] for x in range(0, len(diag_cells), chunk_size)]

                futures = [executor.submit(compute_chunk_tasks,
                       (c, seq1, seq2, match, mismatch, gap, matrix, mode)) for c in chunks]
                
                for f in futures:
                    for i, j, val in f.result():
                        matrix[i, j] = val

    return matrix

#Traceback
def iterative_traceback(matrix, seq1, seq2, match, mismatch, gap, mode):
    """Rebuilds the alignment iteratively"""
    aligned_s1, aligned_s2 = [], []

    if mode == 'global':
        i, j = len(seq1), len(seq2)
    else: 
        i, j = np.unravel_index(np.argmax(matrix), matrix.shape)
    
    score = matrix[i, j]

    while (i > 0 or j > 0):
        if mode == 'local' and matrix[i, j] == 0:
            break

        curr = matrix[i, j]
        #Match/Mismatch

        if i > 0 and j > 0 and curr == matrix[i-1, j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            aligned_s1.append(seq1[i-1])
            aligned_s2.append(seq2[j-1])
            i -= 1; j -= 1
        #Gap in SEQ2
        elif i > 0 and curr == matrix[i-1, j] + gap:
            aligned_s1.append(seq1[i-1])
            aligned_s2.append("-")
            i -= 1
        #Gap in SEQ1 (or edge management j > 0)
        else:
            aligned_s1.append("-")
            aligned_s2.append(seq2[j-1])
            j -= 1

    return "".join(reversed(aligned_s1)), "".join(reversed(aligned_s2)), score

#Execution
if __name__ == '__main__':
    # 1. Load data 
    seq1 = read_fasta("example_data/sequence1.fasta")
    seq2 = read_fasta("example_data/sequence2.fasta")

    #Text example for immediate verification
    #seq1 = "GTCGTAGAATA"
    #seq2 = "CACGTAGCTA"

    if seq1 and seq2:
        match_val, mismatch_pen, gap_pen = 1, -1, -1
        align_mode = 'local' # or 'global'

        print(f"[*] Start calculation {align_mode} su {multiprocessing.cpu_count()} core...")
        start_t = time.time()

        #Step 1: Parallel Matrix Fill
        scoring_matrix = parallel_alignment(seq1, seq2, match_val, mismatch_pen, gap_pen, mode=align_mode)
        
        #Step 2: Traceback (Sequential and fast)
        a1, a2, final_score = iterative_traceback(scoring_matrix, seq1, seq2, match_val, mismatch_pen, gap_pen, mode=align_mode)
        
        duration = time.time() - start_t

        #Output
        print(f"[*] Completed alignment in {duration:.4f} seconds.")
        print(f"[*] Final score: {final_score}")
        
        save_alignment_to_file(a1, a2, final_score, align_mode)
from altools.align.alignment import GlobalAlignment
from altools.dna_sequence import DNASequence

import numpy as np

def align_nw(seq_a: DNASequence, seq_b: DNASequence, match_reward: int = 1, mismatch_penalty: int = -2, gap_penalty: int = -1) -> GlobalAlignment:
    mat = np.zeros(((len(seq_a)+1, len(seq_b)+1)), dtype=int)
    mat[0, :] = [-1*i for i in range(0, len(seq_a)+1)] # seq_a goes along the top
    mat[:, 0] = [-1*i for i in range(0, len(seq_b)+1)] # seq_b goes down the rows

    traceback = {}
    for j in range(1, len(seq_a)+1):
        traceback[(0,j)] = np.array([1])
    for i in range(1, len(seq_b)+1):
        traceback[(i,0)] = np.array([0])

    for i in range(1, len(seq_b)+1):
        for j in range(1, len(seq_a)+1):
            # three options - above value (mat[i-1, j]) + gap penalty
            #                 left value (mat[i, j-1]) + gap penalty
            #                 if they match, diagonal (mat[i-1, j-1]) + match reward
            #                       otherwise, diagonal + mismatch penalty
            top_gap = mat[i-1, j] + gap_penalty
            left_gap = mat[i, j-1] + gap_penalty
            if seq_a[j-1] == seq_b[i-1]:
                matching = mat[i-1, j-1] + match_reward
            else:
                matching = mat[i-1, j-1] + mismatch_penalty
            mat[i,j] = max(top_gap, left_gap, matching)
            scores = np.array([top_gap, left_gap, matching])
            max_indices = np.where(np.max(scores) == scores)
            traceback[(i,j)] = max_indices
    
    start_traceback = (len(seq_b), len(seq_a))


if __name__ == "__main__":
    align_nw(
        DNASequence.from_string("GCATGCG"),
        DNASequence.from_string("GATTACA"),
        match_reward=1,
        mismatch_penalty=-1,
        gap_penalty=-1
    )
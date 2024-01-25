from altools.align.alignment import GlobalAlignment
from altools.dna_sequence import DNASequence

import numpy as np

def align_nw(
    seq_a: DNASequence,
    seq_b: DNASequence,
    match_reward: int = 1,
    mismatch_penalty: int = -2,
    gap_penalty: int = -1,
) -> GlobalAlignment:
    mat = np.zeros(((len(seq_a) + 1, len(seq_b) + 1)), dtype=int)
    mat[0, :] = [-1 * i for i in range(0, len(seq_a) + 1)]  # seq_a goes along the top
    mat[:, 0] = [-1 * i for i in range(0, len(seq_b) + 1)]  # seq_b goes down the rows

    traceback = {(0, 0): np.array([])}
    for j in range(1, len(seq_a) + 1):
        traceback[(0, j)] = np.array([1])
    for i in range(1, len(seq_b) + 1):
        traceback[(i, 0)] = np.array([0])

    for i in range(1, len(seq_b) + 1):
        for j in range(1, len(seq_a) + 1):
            # three options - above value (mat[i-1, j]) + gap penalty
            #                 left value (mat[i, j-1]) + gap penalty
            #                 if they match, diagonal (mat[i-1, j-1]) + match reward
            #                       otherwise, diagonal + mismatch penalty
            top_gap = mat[i - 1, j] + gap_penalty
            left_gap = mat[i, j - 1] + gap_penalty
            if seq_a[j - 1] == seq_b[i - 1]:
                matching = mat[i - 1, j - 1] + match_reward
            else:
                matching = mat[i - 1, j - 1] + mismatch_penalty
            mat[i, j] = max(top_gap, left_gap, matching)
            scores = np.array([top_gap, left_gap, matching])
            max_indices = np.where(np.max(scores) == scores)[0]
            traceback[(i, j)] = max_indices

    traceback_pos = (len(seq_b), len(seq_a))
    # if it is a two, i-1, j-1
    # if it is a 1 i, j-1
    # it it is a 0 i-1, j
    ## what about when there are two options
    ## should we return multiple best alignments? - split cases with multiple into new cases?
    breakpoint()
    traces = [traceback_pos]
    while traceback_pos is not (0, 0):
        pass


def pos_from_traceback(
    traceback: dict[tuple[int, int], np.ndarray], position: tuple[int, int]
) -> list[tuple[int, int]]:
    """
    Computes a new position (or sequence of positions) from a traceback dictionary and a position
    """
    direction_mapping = {0: [-1,0], 1: [0, -1], 2: [-1, -1]}
    positions = []
    for direction_indicator in traceback[position]:
        direction = direction_mapping[direction_indicator]
        new_position = tuple(a + b for (a, b) in zip(position, direction))
        positions.append(new_position)
    return positions

def get_paths(traceback: dict[tuple[int, int], np.ndarray], position: tuple[int, int]) -> list[list[tuple[int, int]]]:
    if position == (0,0):
        return [[]]
    result = []
    for next_path_position in pos_from_traceback(traceback, position):
        sub_paths = get_paths(traceback, next_path_position)
        for sub_path in sub_paths:
            result.append(sub_path + [position])
    return result



if __name__ == "__main__":
    align_nw(
        DNASequence.from_string("GCATGCG"),
        DNASequence.from_string("GATTACA"),
        match_reward=1,
        mismatch_penalty=-1,
        gap_penalty=-1,
    )

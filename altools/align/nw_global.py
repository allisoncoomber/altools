from altools.align.alignment import Alignment
from altools.dna_sequence import DNASequence
from altools.nucleotide import Nucleotide

import numpy as np


def global_align_nw(
    seq_a: DNASequence,
    seq_b: DNASequence,
    match_reward: int = 1,
    mismatch_penalty: int = -2,
    gap_penalty: int = -1,
) -> list[Alignment]:
    mat = np.zeros(((len(seq_b) + 1, len(seq_a) + 1)), dtype=int)
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

    paths = get_paths(traceback, traceback_pos)
    return get_alignments_from_paths(paths, seq_a, seq_b)


def pos_from_traceback(
    traceback: dict[tuple[int, int], np.ndarray], position: tuple[int, int]
) -> list[tuple[int, int]]:
    """
    Computes a new position (or sequence of positions) from a traceback dictionary and a position
    """
    direction_mapping = {0: [-1, 0], 1: [0, -1], 2: [-1, -1]}
    positions = []
    for direction_indicator in traceback[position]:
        direction = direction_mapping[direction_indicator]
        new_position = tuple(a + b for (a, b) in zip(position, direction))
        positions.append(new_position)
    return positions


def get_paths(
    traceback: dict[tuple[int, int], np.ndarray], position: tuple[int, int]
) -> list[list[tuple[int, int]]]:
    if position == (0, 0):
        return [[]]
    result = []
    for next_path_position in pos_from_traceback(traceback, position):
        sub_paths = get_paths(traceback, next_path_position)
        for sub_path in sub_paths:
            result.append(sub_path + [position])
    return result


def get_alignments_from_paths(
    paths: list[list[tuple[int, int]]], seq_a: DNASequence, seq_b: DNASequence
) -> list[Alignment]:
    # start at the first tuple in the path
    # this corresponds to the index+1 of the sequences
    # where i is the seq_b and j is seq_a
    # follow the path, building up both sequences as we go
    # if the indices change the same amount, get the nucleotide at that position from both sequences
    # if there is a left gap, then seq_b has a gap
    # if there is top gap, then seq_a has a gap
    # okay this is three branches, is there a way to make it less branch-y
    # i think yes if instead of making "three cases" we just have one situation
    # where we compare to the previous value each time?
    alignments = []
    for path in paths:
        aligned_a = []
        aligned_b = []
        for idx in range(len(path)):
            i_current, j_current = path[idx]
            if idx == 0:
                difference_in_position = (i_current, j_current)
            else:
                i_previous, j_previous = path[idx - 1]
                difference_in_position = (
                    i_current - i_previous,
                    j_current - j_previous,
                )
            if difference_in_position == (1, 1):
                aligned_a.append(seq_a[j_current - 1])
                aligned_b.append(seq_b[i_current - 1])
            elif difference_in_position == (0, 1):
                aligned_a.append(seq_a[j_current - 1])
                aligned_b.append(Nucleotide.GAP)
            elif difference_in_position == (1, 0):
                aligned_a.append(Nucleotide.GAP)
                aligned_b.append(seq_b[i_current - 1])
        alignments.append(
            Alignment(DNASequence(aligned_a), DNASequence(aligned_b))
        )
    return alignments

import math

import numpy as np
from Bio import SeqIO
from Bio.Align import substitution_matrices


H = -11  # gap opening penalty
G = -1  # gap extension penalty
S = substitution_matrices.load("BLOSUM62")  # scoring matrix


def parse_fasta_file(fasta_file_path: str) -> tuple[np.ndarray, np.ndarray]:
    records = list(SeqIO.parse(fasta_file_path, 'fasta'))
    assert len(records) == 2, "wrong number of records in the provided fasta file"
    return np.array(records[0].seq), np.array(records[1].seq)


def print_matrix(matrix: np.ndarray, seq_a: np.ndarray, seq_b: np.ndarray, pos: list[tuple[int, int]]) -> None:
    height, width, _ = matrix.shape
    for j in range(width - 1):
        for j in range(width):
            print(f" {seq_a[j - 1] if j > 0 else '   '} ", end="                   ")
    print(f"")
    for i in range(height):
        print(f" {seq_b[i - 1] if 0 < i < height - 1 else ' '} ", end="")
        print(f" {seq_b[i - 1] if 0 < i < height else ' '} ", end="")
        for j in range(width):
            if (i, j) in pos:
                print(u"\u001b[1;38;2;255;255;75m", end="")
            print("[", end="")
            for k in range(3):
                print(f"{matrix[i][j][k]:5}", end=" ")
            print("], ", end="")
            if (i, j) in pos:
                print(u"\u001b[0m", end="")
        print("")


def construct_matrix(seq_a: np.ndarray, seq_b: np.ndarray) -> np.array:
    # Initialise Matrix
    m = np.zeros((seq_b.size + 1, seq_a.size + 1, 3))
    m[0][0] = [0, 0, 0]
    for j in range(1, seq_a.size + 1):
        m[0][j] = [-np.inf, (j - 1) * G + H, (j - 1) * G + H]
    for i in range(1, seq_b.size + 1):
        m[i][0] = [(i - 1) * G + H, (i - 1) * G + H, -np.inf]

    # Construct matrix
    for i in range(1, seq_b.size + 1):
        for j in range(1, seq_a.size + 1):
            m[i][j][0] = max(m[i - 1][j][0] + G, m[i - 1][j][1] + H)
            m[i][j][2] = max(m[i][j - 1][2] + G, m[i][j - 1][1] + H)
            m[i][j][1] = max(m[i - 1][j - 1][1] + S[seq_a[j - 1]][seq_b[i - 1]], m[i][j][0], m[i][j][2])

    # Return matrix
    return m


def global_alignment_score(fasta_file_path: str) -> int:
    seq_a, seq_b = parse_fasta_file(fasta_file_path)
    m = construct_matrix(seq_a, seq_b)
    return int(np.amax(m[-1][-1]))


DIRS = [(-1, 0), (-1, -1), (0, -1)]


def global_alignment(fasta_file_path: str) -> tuple[str, str]:
    # Parse sequences from file and construct matrix
    seq_a, seq_b = parse_fasta_file(fasta_file_path)
    m = construct_matrix(seq_a, seq_b)

    # Initialise traceback
    seq_a_aligned = []
    seq_b_aligned = []

    i = m.shape[0] - 1
    j = m.shape[1] - 1
    platform = 1

    # Traceback
    path = []
    while i > 0 or j > 0:
        path.append((i, j))
        print(platform)

        # Decide direction and platform
        i_dirs = platform
        if platform == 1:
            i_dirs = np.argmax([m[i][j][0], m[i - 1][j - 1][1] + S[seq_a[j - 1]][seq_b[i - 1]], m[i][j][2]])

        if i_dirs == 0:
            platform = 0 if m[i - 1][j][0] + G > m[i - 1][j][1] + H else 1
        elif i_dirs == 2:
            platform = 2 if m[i][j - 1][2] + G > m[i][j - 1][1] + H else 1

        seq_a_aligned.append('-' if i_dirs == 0 else seq_a[j - 1])
        seq_b_aligned.append('-' if i_dirs == 2 else seq_b[i - 1])

        i += DIRS[i_dirs][0]
        j += DIRS[i_dirs][1]

    seq_a_aligned = ''.join(np.flip(seq_a_aligned))
    seq_b_aligned = ''.join(np.flip(seq_b_aligned))

    print_matrix(m, seq_a, seq_b, path)

    return seq_a_aligned, seq_b_aligned

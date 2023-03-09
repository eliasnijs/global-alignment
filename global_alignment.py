import numpy as np
from Bio import SeqIO
from Bio.Align import substitution_matrices

H = -11  # Gap opening penalty
G = -1  # Gap extension penalty
S = substitution_matrices.load("BLOSUM62")  # scoring matrix


def print_matrix(matrix: np.ndarray, seq_a: np.ndarray, seq_b: np.ndarray, pos: list[tuple[int, int]]) -> None:
    height, width, _ = matrix.shape
    for j in range(width - 1):
        print(f" {seq_a[j - 1] if j > 0 else '   '} ", end="                   ")
    print(f"")
    for i in range(height):
        print(f" {seq_b[i - 1] if 0 < i < height - 1 else ' '} ", end="")
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


def parse_fasta_file(fasta_file_path: str) -> tuple[np.ndarray, np.ndarray]:
    records = list(SeqIO.parse(fasta_file_path, 'fasta'))
    assert len(records) == 2, "wrong number of records in the provided fasta file"
    return np.array(records[0].seq), np.array(records[1].seq)


def construct_matrix(seq_a: np.ndarray, seq_b: np.ndarray) -> np.array:
    # Initialise Matrix
    m = np.zeros((seq_b.size + 1, seq_a.size + 1, 3))
    m[0][0] = [0, H, H]
    for j in range(1, seq_a.size + 1):
        m[0][j] = [-np.inf, -np.inf, H + j * G]
        m[-1][j] = [-np.inf, -np.inf, -np.inf]
    for i in range(1, seq_b.size + 1):
        m[i][0] = [-np.inf, H + i * G, -np.inf]
        m[i][-1] = [-np.inf, -np.inf, -np.inf]

    # Fill in matrix
    for i in range(1, seq_b.size + 1):
        for j in range(1, seq_a.size + 1):
            i1 = seq_a[j - 1] if j <= seq_a.size else '*'
            i2 = seq_b[i - 1] if i <= seq_b.size else '*'
            m[i][j] = [
                max(m[i - 1][j - 1]) + S[i1][i2],
                max(m[i - 1][j][0] + H, m[i - 1][j][1] + G),
                max(m[i][j - 1][0] + H, m[i][j - 1][2] + G)
            ]

    # Return matrix
    return m


def global_alignment_score(fasta_file_path: str) -> int:
    seq_a, seq_b = parse_fasta_file(fasta_file_path)
    m = construct_matrix(seq_a, seq_b)

    print(m[-2][-2], m[-2][-1])
    print(m[-1][-2], m[-1][-1])

    return int(np.amax(m[-1][-1]))


def global_alignment(fasta_file_path: str) -> tuple[str, str]:
    dirs = np.array([[-1, -1], [0, -1], [-1, 0]])  # Directions in traceback

    seq_a, seq_b = parse_fasta_file(fasta_file_path)
    m = construct_matrix(seq_a, seq_b)

    seq_a_aligned = []
    seq_b_aligned = []

    # Initialise traceback
    i = m.shape[0] - 1
    j = m.shape[1] - 1

    # Traceback
    while i > 1 or j > 1:
        i_dirs = np.argmax([m[i - 1][j - 1][0], m[i - 1][j][2], m[i][j - 1][1]])
        i += dirs[i_dirs][0]
        j += dirs[i_dirs][1]
        seq_a_aligned.append('-' if i_dirs == 2 else seq_a[j - 1])
        seq_b_aligned.append('-' if i_dirs == 1 else seq_b[i - 1])

    # Return aligned sequences
    return ''.join(np.flip(seq_a_aligned)), ''.join(np.flip(seq_b_aligned))
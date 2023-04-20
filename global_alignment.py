import numpy as np
from Bio import SeqIO
from Bio.Align import substitution_matrices

H = np.int16(-11)  # gap opening penalty
G = np.int16(-1)  # gap extension penalty

S = np.array(substitution_matrices.load("BLOSUM62"), dtype=np.int8)  # Scoring matrix

CONV_TABLE = {  # Table for converting letter to index in BLOSUM62
    'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C': 4, 'Q': 5, 'E': 6, 'G': 7, 'H': 8, 'I': 9, 'L': 10, 'K': 11, 'M': 12, 'F': 13,
    'P': 14, 'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19, 'B': 20, 'Z': 21, 'X': 22, '*': 23,
}


def parse_fasta_file(fasta_file_path: str) -> tuple[np.ndarray, np.ndarray]:
    records = list(SeqIO.parse(fasta_file_path, 'fasta'))
    assert len(records) == 2, "wrong number of records in the provided fasta file"
    return np.array(records[0].seq), np.array(records[1].seq)


def sequence_to_indices(seq: np.ndarray):
    index_sequence = np.zeros(len(seq), dtype=np.int8)
    for i, letter in enumerate(seq):
        index_sequence[i] = CONV_TABLE[letter]
    return index_sequence


def construct_matrix(index_seq_a: np.ndarray, index_seq_b: np.ndarray) -> np.array:
    # Initialise matrix
    m = np.zeros((index_seq_b.size + 1, index_seq_a.size + 1, 3), dtype=np.int16)
    m[0][0] = [0, 0, 0]
    for j in range(1, index_seq_a.size + 1):
        m[0][j] = [-32768 / 2, (j - 1) * G + H, (j - 1) * G + H]
    for i in range(1, index_seq_b.size + 1):
        m[i][0] = [(i - 1) * G + H, (i - 1) * G + H, -32768 / 2]

    # Construct matrix
    for i in range(1, index_seq_b.size + 1):
        for j in range(1, index_seq_a.size + 1):
            m[i][j][0] = max(m[i - 1][j][0] + G, m[i - 1][j][1] + H)
            m[i][j][2] = max(m[i][j - 1][2] + G, m[i][j - 1][1] + H)
            m[i][j][1] = max(m[i - 1][j - 1][1] + S[index_seq_a[j - 1]][index_seq_b[i - 1]], m[i][j][0], m[i][j][2])

    # Return matrix
    return m


def global_alignment_score(fasta_file_path: str) -> int:
    seq_a, seq_b = parse_fasta_file(fasta_file_path)
    index_seq_a = sequence_to_indices(seq_a)
    index_seq_b = sequence_to_indices(seq_b)
    m = construct_matrix(index_seq_a, index_seq_b)
    return np.amax(m[-1][-1])


DIRS = [(-1, 0), (-1, -1), (0, -1)]


def global_alignment(fasta_file_path: str) -> tuple[str, str]:
    # Parse sequences from file and construct matrix
    seq_a, seq_b = parse_fasta_file(fasta_file_path)
    index_seq_a = sequence_to_indices(seq_a)
    index_seq_b = sequence_to_indices(seq_b)
    m = construct_matrix(index_seq_a, index_seq_b)

    # Initialise traceback
    seq_a_aligned = []
    seq_b_aligned = []

    i = m.shape[0] - 1
    j = m.shape[1] - 1
    platform = 1

    # Traceback
    while i > 0 or j > 0:
        # Decide direction and platform
        i_dirs = platform
        if platform == 1:
            i_dirs = np.argmax([m[i][j][0], m[i - 1][j - 1][1] + S[index_seq_a[j - 1]][index_seq_b[i - 1]], m[i][j][2]])

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

    return seq_a_aligned, seq_b_aligned

from global_alignment import *

# Helper functions


def print_red():
    # print(u"\u001b[1;38;2;255;75;75m", end="")
    pass


def print_green():
    # print(u"\u001b[1;38;2;75;255;75m", end="")
    pass


def print_normal():
    # print(u"\u001b[0m", end="")
    pass

def test_score(file, correct_score):
    score = global_alignment_score(file)
    if score == correct_score:
        print_green()
        print("> test ok!")
    else:
        print_red()
        print(f"> test failed! Expected {correct_score} but was {score}")
    print_normal()
    print()


def test_sequences(file, correct_sequences):
    seq_a, seq_b = global_alignment(file)

    print()
    print("> Result:")
    if seq_a != correct_sequences[0]:
        print(f"seq_a was incorrect")
        print(f"our seq_a:\t{seq_a}")
        print(f"correct seq_a:\t{correct_sequences[0]}")
    if seq_b != correct_sequences[1]:
        print(f"seq_b was incorrect")
        print(f"our seq_b\t{seq_b}")
        print(f"correct seq_b\t{correct_sequences[1]}")

    if seq_a == correct_sequences[0] and seq_b == correct_sequences[1]:
        print_green()
        print("> test ok!")
    else:
        print_red()
        print("> test failed!")
    print_normal()

# Tests
def test_fasta_01():
    file = "./resources/data01.fasta"

    # Test score
    test_score(file, 8)

    # Test sequences
    correct_sequences = ['PRT---EINS', 'PRTWPSEIN-']
    test_sequences(file, correct_sequences)


def test_fasta_02():
    file = "./resources/data02.faa"

    # Test score
    test_score(file, 144)

    # Test sequences
    correct_sequences = [
        'YHFDVPDCWAHRYWVENPQAIAQME-------QICFNWFPSMMMK-------QPHVFKV---DHHMSCRWLPIRGKKCSSCCTRMRVRTVWE',
        'YHEDV----AHE------DAIAQMVNTFGFVWQICLNQFPSMMMKIYWIAVLSAHVADRKTWSKHMSCRWLPI----ISATCARMRVRTVWE',
    ]
    test_sequences(file, correct_sequences)



























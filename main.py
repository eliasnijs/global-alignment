from global_alignment import *
import tests


def main():
    print(f"////////////////////////////////////////")
    print(f"//// test_fasta_01")
    tests.test_fasta_01()

    print(f"////////////////////////////////////////")
    print(f"//// test_fasta_02")
    tests.test_fasta_02()

    print(f"////////////////////////////////////////")
    print(f"//// test_fasta_06")
    tests.test_fasta_06()

    print(f"////////////////////////////////////////")
    print(f"//// test_fasta_08")
    tests.test_fasta_08()

    print(f"////////////////////////////////////////")
    print(f"//// test_fasta_12")
    tests.test_fasta_12()

    print(f"////////////////////////////////////////")
    print(f"//// test_fasta_17")
    tests.test_fasta_17()


if __name__ == '__main__':
    main()
from global_alignment import *

# Helper functions


def print_red():
    print(u"\u001b[1;38;2;200;75;75m", end="")
    pass


def print_green():
    print(u"\u001b[1;38;2;75;200;75m", end="")
    pass


def print_normal():
    print(u"\u001b[0m", end="")
    pass


def test_score(file, correct_score):
    score = global_alignment_score(file)
    if score == correct_score:
        print_green()
        print("> score test:\t\t[OK]")
    else:
        print(f"Expected {correct_score} but was {score}")
        print_red()
        print(f"> score test:\t\t[failed]")
    print_normal()


def test_sequences(file, correct_sequences):
    seq_a, seq_b = global_alignment(file)
    if seq_a != correct_sequences[0]:
        print_red()
        print(f"seq_a failed:")
        print_normal()
    else:
        print_green()
        print(f"seq_a ok:")
    print_normal()
    print(f"\tour seq_a:\t{seq_a}")
    print(f"\tcorrect seq_a:\t{correct_sequences[0]}")
    if seq_b != correct_sequences[1]:
        print_red()
        print(f"seq_b failed:")
        print_normal()
    else:
        print_green()
        print(f"seq_b ok:")
        print_normal()
    print(f"\tour seq_b\t{seq_b}")
    print(f"\tcorrect seq_b\t{correct_sequences[1]}")
    if seq_a == correct_sequences[0] and seq_b == correct_sequences[1]:
        print_green()
        print("> sequence test\t\t[OK]")
    else:
        print_red()
        print("> sequence test\t\t[failed]")
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
        'YHEDV----AHE------DAIAQMVNTFGFVWQICLNQFPSMMMKIYWIAVLSAHVADRKTWSKHMSCRWLPI----ISATCARMRVRTVWE'
    ]
    test_sequences(file, correct_sequences)


def test_fasta_06():
    file = "./resources/data06.faa"

    correct_sequences = [
        'GKRYCKVMRFDLGIEVNL-TPTVLKHFTMMTHP-NFVTCWF-HKP-QVHE-CHKSVLP-DALWGQMKGNNWATNNIAVNDGFKA'
        '-ETVSFTLCYKTEETTPNMSRLMCLHPLMWEQNEQTAVLVQWPKEWYNVPFEYPM-CLPEDKLCPCQGYQQAEIPMIKWIDE'
        '-FQWWYDMWALHWVSTTGWFNPVNINGLETVCFDMHDMLSQPR'
        '--ADTHADFRGMSQHFPPAHGFRMGNLRVYGHMPGEEGFKYFPHVIKTQARAQMDKMHPNGMLGTIIMPHFDYRCCM-GDTEC-MRMGS',
        '-KRYCKVMRFDL-IEVNLMTPTVLKHTMLMTHPGNTFVCWFDHKQRQVHEDCHKSVLPLDALW-QMKGNNWATNIEAVGDGFKASETVSFTLCVKTE-TTPNM'
        '-RLMCLHPLMWEQNE-TAVL-QWPKEWYHVPFEYQMTCLPEDKLCPCQGYQQAEIPMIKWIDEVFQWWYDMWALHWVSTTGWFNPVNIN'
        '-LPTVCFDMHDMLSEPRFEADTHADGR-MSQHFPPAHGFRMGNLRVYGHMPG--GFK-FPHVIKTQARAQMDIAHPNGMLG-IIMPHFDYRCCMWGEWECHMRMGS'
    ]

    test_score(file, 1184)
    test_sequences(file, correct_sequences)



def test_fasta_08():
    file = "./resources/data08.faa"

    # Test score
    test_score(file, 319)


def test_fasta_12():
    file = "./resources/data12.faa"

    # Test score
    test_score(file, 466)

    # Test sequences
    correct_sequences = [
        'FCMHGYGDVLNACNDIAVCFWINADPVQTTFVIC'
        '-ENQRCKEGKWQNSFKESGHATLFEWWQDMTCQGEYNVMSGPFQNGASIESIKEPLCRMGDQGPTEMTCVLDEVCWGLMSIQATPCSYTPIHIRSM'
        '-EGPVSEAHHVEVWWPNSYLWYYCIRFNADQRCMCN-TKMIQWSKMYEDWSDTKCSQEFSMEWYRKSWF--HAESDK'
        '-LTYVMPMLETDENPWHPMLMPKADAFIRKFQNSNIC-LFPECHLDIAKNYNDKKSPWTGTRKYGTQSLLVTRWHKENVKCSRTKKESCVTHSWSQIAIRPAP'
        '-TQYKLQKAPWPHMSGIDWEAVNIQQRWLIEMVYIYRDHCEHTSWFKHGALDFAF-RCIVWTVDTFPESWC',
        'WMCTSGGGCLNC---IIPCFC--RDPVH----MCPENQRMKLG-WDCSFSLSGHAQLFSLI-DMI-GGECNVMTLSFQNGHASSGSTCMLCKMVDQGTT'
        '-CTMVLVEVCFGDMWIEYP--RVVPQHIASMREIPVFAGWEEVIEYME--LWYIDVGF--EHTCMCNETGQIIWSKMCE-LTCTLC'
        '--EFWMYPYCKPWFESHALHDQNVGYLASMLETDW-PSHPMLEKHFAAY--DFFTVNILEYFKEMHLCIVKNYNIES--WTDTRPYGWQ-VLSTRWMKENC'
        '-CSRPIPWACTTVFWVCIAIRKYPQSQYKQMQKYW--QMGIDCEA-NIQQRPK-EQYYSYRSKCEHTYWNK-SELDMAFHKTIVWTV-IHPESW-'
    ]
    test_sequences(file, correct_sequences)


def test_fasta_17():
    file = "./resources/data17.faa"

    # Test score
    test_score(file, 1695)
















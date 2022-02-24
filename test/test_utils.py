import unittest

from bioseq import utils, DNA, Peptide, Sequence
from test.test_bioseq import TEST_DNA


class TestUtils(unittest.TestCase):
    def test_printAlign(self):
        pass

    def test_parseFasta(self):
        copies = 2
        fasta_text = "\n".join([TEST_DNA] * copies)

        datas = TEST_DNA.splitlines()
        info = datas[0].lstrip(">")
        sequence = "".join(datas[1:]).upper()

        seqs = utils.parseFasta(fasta_text)
        self.assertEqual(len(seqs), copies)

        for seq in seqs:
            self.assertEqual(seq.info, info)
            self.assertEqual(seq, sequence)

    def test_loadFasta(self):
        copies = 2
        with open("test.fasta", "w") as f:
            for _ in range(copies):
                f.write(TEST_DNA)

        datas = TEST_DNA.splitlines()
        info = datas[0].lstrip(">")
        sequence = "".join(datas[1:]).upper()

        seqs = utils.loadFasta("test.fasta")
        self.assertEqual(len(seqs), copies)

        for seq in seqs:
            self.assertEqual(seq.info, info)
            self.assertEqual(seq, sequence)

        import os
        os.remove("test.fasta")

    def test_loadFasta_iter(self):
        copies = 2
        with open("test.fasta", "w") as f:
            for _ in range(copies):
                f.write(TEST_DNA)

        datas = TEST_DNA.splitlines()
        info = datas[0].lstrip(">")
        sequence = "".join(datas[1:]).upper()

        seqs = list(utils.loadFasta("test.fasta", iterator=True))
        self.assertEqual(len(seqs), copies)

        for seq in seqs:
            self.assertEqual(seq.info, info)
            self.assertEqual(seq, sequence)

        import os
        os.remove("test.fasta")

    def test_fetchNCBI(self):
        dna = utils.fetchNCBI("NM_001101.5")
        self.assertEqual(
            dna.info, "NM_001101.5 Homo sapiens actin beta (ACTB), mRNA")
        self.assertIsInstance(dna, DNA)
        self.assertEqual(dna.length, 1812)

        peptide = utils.fetchNCBI("NP_001092.1")
        self.assertEqual(
            peptide.info, "NP_001092.1 actin, cytoplasmic 1 [Homo sapiens]")
        self.assertIsInstance(peptide, Peptide)
        self.assertEqual(peptide.length, 375)

    def test_fetchNCBI_MultiUid(self):
        seq = utils.fetchNCBI(["NM_001101.5", "NP_ENS1111", "NP_001092.1"])
        self.assertEqual(
            seq[0].info, "NM_001101.5 Homo sapiens actin beta (ACTB), mRNA")
        self.assertEqual(
            seq[-1].info, "NP_001092.1 actin, cytoplasmic 1 [Homo sapiens]")
        self.assertIsInstance(seq[0], Sequence)
        self.assertIsInstance(seq[-1], Sequence)

        self.assertEqual(seq[0].length, 1812)
        self.assertEqual(seq[1].length, 375)

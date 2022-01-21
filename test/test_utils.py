import unittest

from bioseq import utils, DNA, Peptide
from test.test_bioseq import TEST_DNA


class TestUtils(unittest.TestCase):
    def test_printAlign(self):
        pass

    def test_loadFasta(self):
        copies = 2
        with open("test.fasta", "w") as f:
            for _ in range(copies):
                f.write(TEST_DNA)

        datas = TEST_DNA.splitlines()
        info = datas[0].lstrip(">")
        sequence = "".join(datas[1:]).upper()

        seqs = list(utils.loadFasta("test.fasta"))
        self.assertEqual(len(seqs), copies)

        for seq in seqs:
            self.assertEqual(seq.info, info)
            self.assertEqual(seq, sequence)

        import os
        os.remove("test.fasta")

    def test_fetchNCBI(self):
        dna = utils.fetchNCBI("NM_001101.5")
        self.assertEqual(dna.info, "NM_001101.5 Homo sapiens actin beta (ACTB), mRNA")
        self.assertIsInstance(dna, DNA)
        self.assertEqual(dna.length, 1812)

        peptide = utils.fetchNCBI("NP_001092.1")
        self.assertEqual(peptide.info, "NP_001092.1 actin, cytoplasmic 1 [Homo sapiens]")
        self.assertIsInstance(peptide, Peptide)
        self.assertEqual(peptide.length, 375)

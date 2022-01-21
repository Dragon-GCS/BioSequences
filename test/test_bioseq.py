from bioseq import DNA, RNA, Peptide
from bioseq._sequence import Sequence
from bioseq.config import AlignmentConfig, setAlignPara, MW
import unittest


# NM_001101.5, 85..1212
TEST_DNA = """\
>Homo sapiens actin beta (ACTB), mRNA
atggatgatgatatcgccgcgctcgtcgtcgacaacggctccggcatgtgcaaggccggcttcgcgggcg
acgatgccccccgggccgtcttcccctccatcgtggggcgccccaggcaccagggcgtgatggtgggcat
gggtcagaaggattcctatgtgggcgacgaggcccagagcaagagaggcatcctcaccctgaagtacccc
atcgagcacggcatcgtcaccaactgggacgacatggagaaaatctggcaccacaccttctacaatgagc
tgcgtgtggctcccgaggagcaccccgtgctgctgaccgaggcccccctgaaccccaaggccaaccgcga
gaagatgacccagatcatgtttgagaccttcaacaccccagccatgtacgttgctatccaggctgtgcta
tccctgtacgcctctggccgtaccactggcatcgtgatggactccggtgacggggtcacccacactgtgc
ccatctacgaggggtatgccctcccccatgccatcctgcgtctggacctggctggccgggacctgactga
ctacctcatgaagatcctcaccgagcgcggctacagcttcaccaccacggccgagcgggaaatcgtgcgt
gacattaaggagaagctgtgctacgtcgccctggacttcgagcaagagatggccacggctgcttccagct
cctccctggagaagagctacgagctgcctgacggccaggtcatcaccattggcaatgagcggttccgctg
ccctgaggcactcttccagccttccttcctgggcatggagtcctgtggcatccacgaaactaccttcaac
tccatcatgaagtgtgacgtggacatccgcaaagacctgtacgccaacacagtgctgtctggcggcacca
ccatgtaccctggcattgccgacaggatgcagaaggagatcactgccctggcacccagcacaatgaagat
caagatcattgctcctcctgagcgcaagtactccgtgtggatcggcggctccatcctggcctcgctgtcc
accttccagcagatgtggatcagcaagcaggagtatgacgagtccggcccctccatcgtccaccgcaaat
gcttctag
"""

TEST_PEPTIDE = \
"MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTL" \
"KYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTEAPLNPKANREKMTQIMFETFNTPAMYVAIQ" \
"AVLSLYASGRTTGIVMDSGDGVTHTVPIYEGYALPHAILRLDLAGRDLTDYLMKILTERGYSFTTTAERE" \
"IVRDIKEKLCYVALDFEQEMATAASSSSLEKSYELPDGQVITIGNERFRCPEALFQPSFLGMESCGIHET" \
"TFNSIMKCDVDIRKDLYANTVLSGGTTMYPGIADRMQKEITALAPSTMKIKIIAPPERKYSVWIGGSILA" \
"SLSTFQQMWISKQEYDESGPSIVHRKCF"

class TestBioseq(unittest.TestCase):
    def test_seq_add(self):
        self.assertEqual(DNA("ATCG"), DNA("ATCG"))
        self.assertEqual(DNA("ATCG") + DNA("ATCG"), DNA("ATCGATCG"))
        with self.assertRaises(TypeError):
            DNA("ATCG") + RNA("AUCG")

    def test_seq_align(self):
        seq_a = Sequence("ATCG")
        seq_b = Sequence("ATCGATCG")
        target_score = 4 * AlignmentConfig.MATCH + \
                       3 * AlignmentConfig.GAP_EXTEND + \
                       AlignmentConfig.GAP_OPEN

        self.assertEqual(seq_a.align(seq_b), 
                        ("ATCG----", "ATCGATCG", target_score))

    def test_seq_find(self):
        self.assertEqual(Sequence("ATCGATCG").find("CGAT")[0], 2)

    def test_seq_mutation(self):
        seq = Sequence("ATCG")
        self.assertEqual(seq.mutation([0, 1, 2], "A"), "AAAG")
        self.assertEqual(seq.mutation("A", "C"), "CCCG")
        self.assertEqual(seq.mutation(0, "ATC"), "ATCG")

    def test_seq_trans(self):
        seq = Sequence()
        self.assertIsInstance(seq.toDNA(), DNA)
        self.assertIsInstance(seq.toRNA(), RNA)
        self.assertIsInstance(seq.toPeptide(), Peptide)


class TestDNA(unittest.TestCase):
    def setUp(self):
        self.dna = DNA("".join(TEST_DNA.splitlines()[1:]))

    def test_GC(self):
        self.assertEqual(DNA("ATCG").GC, 0.5)

    def test_weight(self):
        self.assertEqual(RNA("AUCG").weight, sum(MW["RNA_MW"].values()) + 18)
        self.assertEqual(DNA("ATCG").weight, sum(MW["DNA_MW"].values()) + 18)


    def test_compostion(self):
        compostion = DNA("ATTCCCGGGGG").composition
        self.assertEqual(compostion["A"], 1)
        self.assertEqual(compostion["C"], 3)
        self.assertEqual(compostion["G"], 5)
        self.assertEqual(compostion["T"], 2)

    def test_reversed(self):
        seq = DNA("ATCG")
        self.assertEqual(seq.reversed, "GCTA")
        seq.reverse()
        self.assertEqual(seq, "GCTA")

    def test_complemented(self):
        seq = DNA("ATCG")
        self.assertEqual(seq.complemented, "CGAT")
        seq.complement()
        self.assertEqual(seq, "CGAT")

    def test_traslate(self):
        self.assertEqual(str(self.dna.translate()), str(self.dna).replace("T", "U"))

    def test_getOrf(self):
        self.assertTrue(not self.dna.orf)
        self.assertEqual(self.dna.getOrf()[0], self.dna.translate())
        self.assertEqual(self.dna.orf[0], self.dna.translate())

    def test_transcript(self):
        self.assertTrue(not self.dna.peptide)
        self.assertEqual(self.dna.transcript()[0], TEST_PEPTIDE)
        self.assertEqual(self.dna.peptide[0], TEST_PEPTIDE)

    def test_reset(self):
        attrs = [attr for attr in self.dna.__dict__ 
                if attr not in ["_seq", "info"]]
        for attr in attrs:
            setattr(self.dna, attr, "test")
        self.dna.reset_cache()
        for attr in attrs:
            self.assertTrue(not getattr(self.dna, attr))


class TestPeptide(unittest.TestCase):
    def setUp(self) -> None:
        self.pep = Peptide(TEST_PEPTIDE)

    def test_pI(self):
        pI = self.pep.pI
        self.assertGreaterEqual(self.pep.pI, 0)
        self.assertLessEqual(self.pep.pI, 7)
        self.assertTrue(self.pep._pI)

    def test_Hphob(self):
        import sys
        win_size = 6
        self.assertEqual(len(self.pep.getHphob(window_size=win_size)), self.pep.length - win_size)
        self.assertTrue(self.pep._Hphob_list)

    def test_reset(self):
        attrs = [attr for attr in self.pep.__dict__ 
                if attr not in ["_seq", "info"]]
        for attr in attrs:
            setattr(self.pep, attr, "test")
        self.pep.reset_cache()
        for attr in attrs:
            self.assertTrue(not getattr(self.pep, attr))
    
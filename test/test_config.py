import unittest

from bioseq import Sequence, config
from bioseq.config import AlignmentConfig, setAlignPara, setStartCoden


class TestUtils(unittest.TestCase):
    def test_alignConfig(self):
        seq_a = Sequence("ATCG")
        seq_b = Sequence("ATCGATCG")

        setAlignPara(5, -4, -4, -4)
        target_score = 4 * AlignmentConfig.MATCH + \
                       3 * AlignmentConfig.GAP_EXTEND + \
                       AlignmentConfig.GAP_OPEN
        self.assertEqual(target_score, 4)
        self.assertEqual(seq_a.align(seq_b), 
                        ("ATCG----", "ATCGATCG", target_score))

    def test_setStartCoden(self):
        coden = ["AAA", "TTT"]
        default = config.START_CODON
        setStartCoden(coden)
        self.assertListEqual(coden, config.START_CODON)
        setStartCoden()
        self.assertListEqual(default, config.START_CODON)

        with self.assertRaises(ValueError):
            # length of each start coden must equals to three
            setStartCoden("1234")

        with self.assertRaises(ValueError):
            setStartCoden(1234)
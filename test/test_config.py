import unittest

from bioseq import Sequence, config
from bioseq.config import AlignmentConfig


class TestUtils(unittest.TestCase):
    def test_alignConfig(self):
        seq_a = Sequence("ATCG")
        seq_b = Sequence("ATCGATCG")

        AlignmentConfig.MATCH,\
        AlignmentConfig.MISMATCH, \
        AlignmentConfig.GAP_OPEN, \
        AlignmentConfig.GAP_EXTEND = (5, -4, -4, -4)

        target_score = 4 * AlignmentConfig.MATCH + \
                       3 * AlignmentConfig.GAP_EXTEND + \
                       AlignmentConfig.GAP_OPEN
        self.assertEqual(target_score, 4)
        self.assertEqual(seq_a.align(seq_b), 
                        ("ATCG----", "ATCGATCG", target_score))
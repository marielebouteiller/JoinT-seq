#!/usr/bin/env python3

import unittest

# Some naive algorithms to help align sequences or count mismatches


def approximate_match(seqA: str, seqB: str, max_mismatchs: int) -> bool:
    """ Match two strings, allowing up to n mismatches.
    This only tests from the starting position
    """

    n_mismatch = 0

    # Go through all characters one by one
    for a, b in zip(seqA, seqB):
        if a != b:
            n_mismatch += 1

        # Early exit if there are already too many mismatches
        if n_mismatch > max_mismatchs:
            return False

    return True


class TestApproxMatch(unittest.TestCase):
    def test_perfect_match(self):
        self.assertTrue(approximate_match("aaaaa", "aaaa", 0))

    def test_one_mismatch(self):
        self.assertTrue(approximate_match("aataa", "aaaaa", 1))
        self.assertFalse(approximate_match("aatat", "aaaaa", 1))

    def test_different_lengths(self):
        self.assertTrue(approximate_match("atcg", "atc", 0))
        self.assertTrue(approximate_match("a", "gcatc", 1))
        self.assertFalse(approximate_match("act", "gcatc", 1))


if __name__ == "__main__":
    unittest.main()

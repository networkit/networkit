#!/usr/bin/env python3
import math
import unittest

from networkit.helpers import ranked


class TestHelpers(unittest.TestCase):
	def setUp(self):
		self.nan = float("nan")

	def testRankedEmptySample(self):
		self.assertEqual(ranked([]), [])

	def testRankedSample(self):
		self.assertEqual(ranked([3.0, 1.0, 2.0, 1.0]), [4.0, 1.5, 3.0, 1.5])

	def testRankedNanOnly(self):
		"""A single NaN must not crash; its rank is NaN."""
		ranks = ranked([self.nan])
		self.assertEqual(len(ranks), 1)
		self.assertTrue(math.isnan(ranks[0]))

	def testRankedNanWithFiniteValues(self):
		"""NaN entries get NaN ranks; finite values keep average ranks."""
		ranks = ranked([1.0, self.nan, 2.0])
		self.assertTrue(math.isnan(ranks[1]))
		self.assertEqual(ranks[0], 1.0)
		self.assertEqual(ranks[2], 2.0)

	def testRankedNanOrderIndependent(self):
		"""Ranks must follow their elements under permutation: ties share
		average ranks, NaN positions get NaN ranks, infinities rank at the
		extremes, and each element keeps its own rank regardless of order."""
		base = [2.0, self.nan, 1.0, 2.0, float("-inf"), float("inf")]
		expected_by_value = {
			float("-inf"): 1.0,
			1.0: 2.0,
			2.0: 3.5,  # tie: average of ranks 3 and 4
			float("inf"): 5.0,
		}
		permutations = [
			base,
			[float("inf"), 2.0, self.nan, 1.0, float("-inf"), 2.0],
			[2.0, 2.0, 1.0, float("inf"), float("-inf"), self.nan],
		]

		reference = None
		for sample in permutations:
			ranks = ranked(sample)
			self.assertEqual(len(ranks), len(sample))
			# per-element correspondence: every rank must match the value
			# actually sitting at that position
			for value, rank in zip(sample, ranks):
				if value != value:  # NaN
					self.assertTrue(math.isnan(rank))
				else:
					self.assertEqual(rank, expected_by_value[value])
			# elementwise permutation invariant: (value, rank) pairs, with
			# NaN ranks excluded, must not depend on input order
			pairs = sorted(
				(value, rank)
				for value, rank in zip(sample, ranks)
				if not math.isnan(rank)
			)
			if reference is None:
				reference = pairs
			else:
				self.assertEqual(pairs, reference)

	def testRankedAllNan(self):
		ranks = ranked([self.nan, self.nan])
		self.assertEqual(len(ranks), 2)
		self.assertTrue(all(math.isnan(r) for r in ranks))

	def testRankedFiniteUnchanged(self):
		"""The existing average-rank behavior for finite input is preserved."""
		self.assertEqual(ranked([1.0, 2.0, 2.0]), [1.0, 2.5, 2.5])
		self.assertEqual(ranked([3.0, 1.0, 2.0]), [3.0, 1.0, 2.0])


if __name__ == "__main__":
	unittest.main()

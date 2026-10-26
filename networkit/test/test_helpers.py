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
		"""The same multiset must give the same ranks regardless of input order."""
		sample_a = [self.nan, 1.0, 2.0]
		sample_b = [2.0, self.nan, 1.0]
		ranks_a = ranked(sample_a)
		ranks_b = ranked(sample_b)

		def summary(ranks):
			return (
				sorted(r for r in ranks if not math.isnan(r)),
				sum(1 for r in ranks if math.isnan(r)),
			)

		self.assertEqual(summary(ranks_a), summary(ranks_b))
		self.assertEqual(summary(ranks_a), ([1.0, 2.0], 1))

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

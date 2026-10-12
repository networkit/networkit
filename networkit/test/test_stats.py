#!/usr/bin/env python3
import unittest

import networkit as nk


class TestStats(unittest.TestCase):
	def testGiniTypical(self):
		self.assertAlmostEqual(nk.stats.gini([1, 2, 3]), 0.2222222222222222)

	def testGiniEqualValues(self):
		self.assertEqual(nk.stats.gini([5, 5, 5]), 0.0)

	def testGiniAllZero(self):
		self.assertEqual(nk.stats.gini([0, 0, 0]), 0.0)

	def testGiniEmpty(self):
		self.assertEqual(nk.stats.gini([]), 0.0)


if __name__ == "__main__":
	unittest.main()

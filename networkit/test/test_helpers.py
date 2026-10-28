#!/usr/bin/env python3
import unittest

from networkit.helpers import ranked


class TestHelpers(unittest.TestCase):
	def testRankedEmptySample(self):
		self.assertEqual(ranked([]), [])

	def testRankedSample(self):
		self.assertEqual(ranked([3.0, 1.0, 2.0, 1.0]), [4.0, 1.5, 3.0, 1.5])


if __name__ == "__main__":
	unittest.main()

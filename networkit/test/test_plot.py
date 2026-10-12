#!/usr/bin/env python3
import unittest
from unittest import mock

import networkit as nk


class TestPlot(unittest.TestCase):
	def testConnectedComponentsSizesEmptyGraphRaises(self):
		"""Graph(0) has no components; reject before max() on empty data."""
		with mock.patch.object(nk.plot, "have_seaborn", True), mock.patch.object(
			nk.plot, "have_plt", True
		), mock.patch.object(nk.plot, "seaborn", create=True), mock.patch.object(
			nk.plot, "plt", create=True
		):
			with self.assertRaises(ValueError) as ctx:
				nk.plot.connectedComponentsSizes(nk.Graph(0))
		self.assertIn("no components", str(ctx.exception))


if __name__ == "__main__":
	unittest.main()

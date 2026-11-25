#!/usr/bin/env python3
import unittest
from unittest import mock

import networkit as nk
import networkit.plot


class TestPlot(unittest.TestCase):
	def mockPlotting(self):
		return mock.patch.multiple(
			nk.plot,
			create=True,
			have_seaborn=True,
			have_plt=True,
			seaborn=mock.DEFAULT,
			plt=mock.DEFAULT,
		)

	def testConnectedComponentsSizesPlotsComponents(self):
		G = nk.Graph(5)
		G.addEdge(0, 1)
		G.addEdge(1, 2)
		G.addEdge(3, 4)
		for relativeSizes in (True, False):
			with self.subTest(relativeSizes=relativeSizes), self.mockPlotting() as mocks:
				nk.plot.connectedComponentsSizes(G, relativeSizes)
				mocks["plt"].pie.assert_called_once()
				args, kwargs = mocks["plt"].pie.call_args
				self.assertEqual(args[0], [3, 2])
				self.assertEqual(kwargs["explode"], [0.1, 0])
				if relativeSizes:
					self.assertEqual(kwargs["autopct"](60), "60%")
					self.assertEqual(kwargs["autopct"](5), "")
				else:
					self.assertEqual(kwargs["autopct"](60), "3")

	def testConnectedComponentsSizesEmptyGraphRaises(self):
		"""Graph(0) has no components; reject before max() on empty data."""
		with self.mockPlotting(), self.assertRaises(ValueError) as ctx:
			nk.plot.connectedComponentsSizes(nk.Graph(0))
		self.assertIn("no components", str(ctx.exception))

	def testConnectedComponentsSizesAllNodesDeletedRaises(self):
		"""Graph with all nodes deleted (n=0, upperNodeIdBound=3) has no components."""
		G = nk.Graph(3)
		for u in range(3):
			G.removeNode(u)
		self.assertEqual(G.numberOfNodes(), 0)
		self.assertEqual(G.upperNodeIdBound(), 3)
		with self.mockPlotting(), self.assertRaises(ValueError) as ctx:
			nk.plot.connectedComponentsSizes(G)
		self.assertIn("no components", str(ctx.exception))


if __name__ == "__main__":
	unittest.main()

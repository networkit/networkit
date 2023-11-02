#!/usr/bin/env python3

import unittest

import networkit as nk

def binom3(n):
	return (n*(n-1)*(n-2)) // 6

def binom4(n):
	return (binom3(n) * (n-3)) // 4

class Test_GraphletsCounter(unittest.TestCase):
	def setUp(self):
		V = self.V = 5
		# C_5 \boxtimes C_5
		self.G = nk.Graph(V*V, False, False)
		for u in range(V):
			v = (u+1) % V
			for w in range(V):
				x = (w+1) % V
				uw = u*V + w
				ux = u*V + x
				vw = v*V + w
				vx = v*V + x
				self.G.addEdge(uw, ux)
				self.G.addEdge(uw, vw)
				self.G.addEdge(uw, vx)
				self.G.addEdge(vw, ux)

	def test_3graphlets(self):
		counter = nk.subgraphs.GraphletsCounter(self.G, 3)
		counts = counter.run().getCounts()
		self.assertIsInstance(counts, list)
		Vsquared = self.V*self.V
		expected = [
			4*Vsquared,
			16*Vsquared,
			4*Vsquared*(Vsquared-13),
			binom3(Vsquared) - 4*Vsquared*(Vsquared-8)
		]
		self.assertEqual(counts, expected)

	def test_4graphlets(self):
		counter = nk.subgraphs.GraphletsCounter(self.G, 4)
		counts = counter.run().getCounts()
		self.assertIsInstance(counts, list)
		Vsquared = self.V*self.V
		expected = [
			4*Vsquared,
			16*Vsquared,
			4*Vsquared*(Vsquared-13),
			binom3(Vsquared) - 4*Vsquared*(Vsquared-8),
			Vsquared,
			8*Vsquared,
			28*Vsquared,
			Vsquared,
			8*Vsquared,
			64*Vsquared,
			4*Vsquared*(Vsquared-15),
			4*Vsquared*(4*Vsquared-69),
			Vsquared*(8*Vsquared-143),
			2*Vsquared*(Vsquared**2 - 35*Vsquared + 326),
		]
		expected.append(binom3(Vsquared) + binom4(Vsquared) - sum(expected))
		self.assertEqual(counts, expected)

	def test_exception(self):
		counter = nk.subgraphs.GraphletsCounter(self.G, 0)
		self.assertRaises(RuntimeError, counter.run)

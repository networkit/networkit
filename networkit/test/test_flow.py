#!/usr/bin/env python3
import unittest

import networkit as nk

class TestFlow(unittest.TestCase):

	def setUp(self):
		self.L = nk.readGraph("input/looptest1.gml", nk.Format.GML) #without self-loops

	def testEdmondsKarp(self):
		G = nk.Graph(5)
		G.addEdge(0,1)
		G.addEdge(1,2)
		G.addEdge(1,3)
		G.addEdge(3,4)
		G.addEdge(2,3)
		G.indexEdges()
		EKL = nk.flow.EdmondsKarp(G, 0, 1)
		EKL.run()
		self.assertEqual(EKL.getFlow(0, 1), 1.0)
		self.assertEqual(EKL.getFlow(1), 0.0)
		self.assertEqual(EKL.getMaxFlow(), 1.0)
		self.assertListEqual(EKL.getSourceSet(), [0])
		self.assertListEqual(EKL.getFlowVector(), [-1.0, 0.0, 0.0, 0.0, 0.0])

	def testDinicConstructorPreconditions(self):
		# Undirected graph
		G_undirected = nk.Graph(2, weighted=True, directed=False)
		with self.assertRaises(RuntimeError) as ctx1:
			nk.flow.Dinic(G_undirected, 0, 1)
		self.assertEqual(str(ctx1.exception), "Dinic algorithm requires directed graph!")

		# Unweighted graph
		G_unweighted = nk.Graph(2, weighted=False, directed=True)
		with self.assertRaises(RuntimeError) as ctx2:
			nk.flow.Dinic(G_unweighted, 0, 1)
		self.assertEqual(str(ctx2.exception), "Dinic algorithm requires weighted graph!")

		# Same source and target
		G_same = nk.Graph(2, weighted=True, directed=True)
		with self.assertRaises(RuntimeError) as ctx3:
			nk.flow.Dinic(G_same, 0, 0)
		self.assertEqual(
			str(ctx3.exception),
			"Dinic algorithm requires `source` and `target` node to be different!",
		)

	def testDinicRunNotCalledThrows(self):
		G = nk.Graph(2, weighted=True, directed=True)
		algo = nk.flow.Dinic(G, 0, 1)
		with self.assertRaises(RuntimeError) as ctx:
			algo.getMaxFlow()
		self.assertEqual(str(ctx.exception), "Error, run must be called first")

	def testDinicNegativeWeightsThrows(self):
		G = nk.Graph(3, weighted=True, directed=True)
		G.addEdge(0, 1, 1.0)
		G.addEdge(1, 2, -1.0)
		algo = nk.flow.Dinic(G, 0, 2)
		with self.assertRaises(RuntimeError):
			algo.run()

	def testDinicThreeDisjointPathsWithParallelEdges(self):
		"""
		Combines:
		- multiple edge-disjoint s–t paths
		- parallel edges + indexEdges()
		Expected max flow: 3 (three unit-capacity branches)
		"""
		G = nk.Graph(5, weighted=True, directed=True)
		# branch via 1
		G.addEdge(0, 1, 1.0)
		G.addEdge(1, 0, 1.0)  # back edge
		G.addEdge(1, 4, 1.0)
		# branch via 2
		G.addEdge(0, 2, 1.0)
		G.addEdge(2, 0, 1.0)  # back edge
		G.addEdge(2, 4, 1.0)
		# branch via 3
		G.addEdge(0, 3, 1.0)
		G.addEdge(3, 0, 1.0)  # back edge
		G.addEdge(3, 4, 1.0)
		G.indexEdges()

		algo = nk.flow.Dinic(G, 0, 4)
		algo.run()
		self.assertAlmostEqual(algo.getMaxFlow(), 3.0, places=12)

	def testDinicDisconnectedGraphs(self):
		"""
		Source and sink in different connected components → max flow must be zero.
		"""
		G = nk.Graph(7, weighted=True, directed=True)
		# component 1
		G.addEdge(0, 1, 10.0)
		G.addEdge(1, 2, 5.0)
		G.addEdge(2, 3, 7.0)
		# component 2
		G.addEdge(4, 5, 11.0)
		G.addEdge(5, 6, 10.0)

		algo = nk.flow.Dinic(G, 0, 5)
		algo.run()
		self.assertAlmostEqual(algo.getMaxFlow(), 0.0, places=12)

	def testDinicNumericalStabilityDecimalSplits(self):
		"""
		Parallel s–t paths with small fractional capacities
		plus a sub-epsilon direct edge that should be ignored
		by the tolerance gating.
		Total incoming capacity at sink should be 1.0.
		"""
		G = nk.Graph(7, weighted=True, directed=True)
		G.addEdge(0, 1, 1.0)
		G.addEdge(1, 2, 0.1)
		G.addEdge(2, 6, 0.1)
		G.addEdge(1, 3, 0.2)
		G.addEdge(3, 6, 0.2)
		G.addEdge(1, 4, 0.3)
		G.addEdge(4, 6, 0.3)
		G.addEdge(1, 5, 0.4)
		G.addEdge(5, 6, 0.4)

		# sub-epsilon direct edge that must effectively be ignored
		G.addEdge(0, 6, 1e-18)

		algo = nk.flow.Dinic(G, 0, 6)
		algo.run()
		self.assertAlmostEqual(algo.getMaxFlow(), 1.0, places=12)

class TestSuccessiveShortestPathMinCostFlow(unittest.TestCase):
	CAPACITY = "capacity"
	SUPPLY = "supply"

	def _check_flow_conservation(self, G, flow, supply_attr):
		"""
		For every node u:
		  outgoing_flow(u) - incoming_flow(u) == supply[u]
		"""
		imbalance = [0.0] * G.numberOfNodes()

		# Iterate edges; Python Graph exposes iterEdges() yielding (u,v) for directed graphs.
		# If your NetworKit version differs, replace with the appropriate edge iterator.
		for (u, v) in G.iterEdges():
			eid = G.edgeId(u, v)
			f = flow[eid]
			imbalance[u] += f
			imbalance[v] -= f

		for u in range(G.numberOfNodes()):
			self.assertAlmostEqual(imbalance[u], supply_attr[u], places=12)

	def _make_graph(self, n, *, weighted=True, directed=True, index_edges=True):
		G = nk.Graph(n, weighted=weighted, directed=directed)
		if index_edges:
			G.indexEdges()
		return G

	def _add_edge_cost_capacity(self, G, capacities, u, v, cost, cap):
		G.addEdge(u, v, cost)
		eid = G.edgeId(u, v)
		capacities[eid] = cap
		return eid

	# ---- constructor / precondition checks ----

	def testConstructorThrowsForUndirected(self):
		G = self._make_graph(2, weighted=True, directed=False, index_edges=True)
		with self.assertRaises(RuntimeError) as ctx:
			nk.flow.SuccessiveShortestPathMinCostFlow(G, self.CAPACITY, self.SUPPLY)
		self.assertEqual(str(ctx.exception), "SuccessiveShortestPathMinCostFlow: Graph must be directed.")

	def testConstructorThrowsForUnweighted(self):
		G = self._make_graph(2, weighted=False, directed=True, index_edges=True)
		with self.assertRaises(RuntimeError) as ctx:
			nk.flow.SuccessiveShortestPathMinCostFlow(G, self.CAPACITY, self.SUPPLY)
		self.assertEqual(str(ctx.exception), "SuccessiveShortestPathMinCostFlow: Graph must be weighted.")

	def testConstructorThrowsForUnindexedEdges(self):
		G = self._make_graph(2, weighted=True, directed=True, index_edges=False)
		with self.assertRaises(RuntimeError) as ctx:
			nk.flow.SuccessiveShortestPathMinCostFlow(G, self.CAPACITY, self.SUPPLY)
		self.assertEqual(str(ctx.exception), "SuccessiveShortestPathMinCostFlow: Graph edges must be indexed.")

	def testConstructorThrowsForMissingEdgeAttribute(self):
		G = self._make_graph(2, weighted=True, directed=True, index_edges=True)
		G.attachNodeAttribute(self.SUPPLY, float)

		with self.assertRaises(RuntimeError) as ctx:
			nk.flow.SuccessiveShortestPathMinCostFlow(G, self.CAPACITY, self.SUPPLY)
		self.assertEqual(
			str(ctx.exception),
			f"SuccessiveShortestPathMinCostFlow: Provided edge attribute '{self.CAPACITY}' not found.",
		)

	def testConstructorThrowsForMissingNodeAttribute(self):
		G = self._make_graph(2, weighted=True, directed=True, index_edges=True)
		G.attachEdgeAttribute(self.CAPACITY, float)

		with self.assertRaises(RuntimeError) as ctx:
			nk.flow.SuccessiveShortestPathMinCostFlow(G, self.CAPACITY, self.SUPPLY)
		self.assertEqual(
			str(ctx.exception),
			f"SuccessiveShortestPathMinCostFlow: Provided node attribute '{self.SUPPLY}' not found.",
		)

	def testConstructorThrowsForNegativeCapacity(self):
		G = self._make_graph(2, weighted=True, directed=True, index_edges=True)
		capacities = G.attachEdgeAttribute(self.CAPACITY, float)
		supply = G.attachNodeAttribute(self.SUPPLY, float)

		G.addEdge(0, 1, 1.0)
		eid = G.edgeId(0, 1)
		capacities[eid] = -5.0
		supply[0] = 0.0
		supply[1] = 0.0

		with self.assertRaises(RuntimeError) as ctx:
			nk.flow.SuccessiveShortestPathMinCostFlow(G, self.CAPACITY, self.SUPPLY)
		self.assertEqual(str(ctx.exception), "SuccessiveShortestPathMinCostFlow: Capacities must be non-negative.")

	def testConstructorThrowsForUnbalancedSupplies(self):
		G = self._make_graph(2, weighted=True, directed=True, index_edges=True)
		capacities = G.attachEdgeAttribute(self.CAPACITY, float)
		supply = G.attachNodeAttribute(self.SUPPLY, float)

		G.addEdge(0, 1, 1.0)
		eid = G.edgeId(0, 1)
		capacities[eid] = 5.0
		supply[0] = -5.0
		supply[1] = 4.0

		with self.assertRaises(RuntimeError) as ctx:
			nk.flow.SuccessiveShortestPathMinCostFlow(G, self.CAPACITY, self.SUPPLY)
		self.assertEqual(
			str(ctx.exception),
			"SuccessiveShortestPathMinCostFlow: Sum of node supplies and demands does not add up to zero.",
		)

	# ---- “run not called” checks ----

	def testRunNotCalledGetTotalCostThrows(self):
		G = self._make_graph(3, weighted=True, directed=True, index_edges=True)
		G.attachEdgeAttribute(self.CAPACITY, float)
		supply = G.attachNodeAttribute(self.SUPPLY, float)
		for u in range(3):
			supply[u] = 0.0

		solver = nk.flow.SuccessiveShortestPathMinCostFlow(G, self.CAPACITY, self.SUPPLY)
		with self.assertRaises(RuntimeError) as ctx:
			solver.getTotalCost()
		self.assertEqual(str(ctx.exception), "Error, run must be called first")

	def testRunNotCalledGetFlowThrows(self):
		G = self._make_graph(3, weighted=True, directed=True, index_edges=True)
		G.attachEdgeAttribute(self.CAPACITY, float)
		supply = G.attachNodeAttribute(self.SUPPLY, float)
		for u in range(3):
			supply[u] = 0.0

		solver = nk.flow.SuccessiveShortestPathMinCostFlow(G, self.CAPACITY, self.SUPPLY)
		with self.assertRaises(RuntimeError) as ctx:
			solver.getFlow()
		self.assertEqual(str(ctx.exception), "Error, run must be called first")

	# ---- corner cases ----

	def testRunThrowsOnNegativeCostCycle(self):
		G = self._make_graph(3, weighted=True, directed=True, index_edges=True)
		capacities = G.attachEdgeAttribute(self.CAPACITY, float)
		supply = G.attachNodeAttribute(self.SUPPLY, float)
		for u in range(3):
			supply[u] = 0.0

		self._add_edge_cost_capacity(G, capacities, 0, 1,  1.0, 3.0)
		self._add_edge_cost_capacity(G, capacities, 1, 2,  1.0, 3.0)
		self._add_edge_cost_capacity(G, capacities, 2, 0, -5.0, 3.0)

		solver = nk.flow.SuccessiveShortestPathMinCostFlow(G, self.CAPACITY, self.SUPPLY)
		with self.assertRaises(RuntimeError) as ctx:
			solver.run()
		self.assertEqual(
			str(ctx.exception),
			"SuccessiveShortestPathMinCostFlow: negative-cost cycle in residual graph",
		)

	def testRunThrowsWhenSuppliesCannotBeSatisfied(self):
		G = self._make_graph(3, weighted=True, directed=True, index_edges=True)
		capacities = G.attachEdgeAttribute(self.CAPACITY, float)
		supply = G.attachNodeAttribute(self.SUPPLY, float)

		supply[0] = 10.0
		supply[1] = -5.0
		supply[2] = -5.0

		self._add_edge_cost_capacity(G, capacities, 0, 1, 1.0, 3.0)
		self._add_edge_cost_capacity(G, capacities, 1, 2, 1.0, 3.0)
		self._add_edge_cost_capacity(G, capacities, 2, 0, 2.0, 3.0)

		solver = nk.flow.SuccessiveShortestPathMinCostFlow(G, self.CAPACITY, self.SUPPLY)
		with self.assertRaises(RuntimeError) as ctx:
			solver.run()
		self.assertEqual(
			str(ctx.exception),
			"SuccessiveShortestPathMinCostFlow: unable to satisfy all supplies/demands",
		)

	# ---- functional tests ----
	def testSimple5NodeGraph(self):
		G = self._make_graph(5, weighted=True, directed=True, index_edges=True)
		capacities = G.attachEdgeAttribute(self.CAPACITY, float)
		supply = G.attachNodeAttribute(self.SUPPLY, float)

		for u in range(5):
			supply[u] = 0.0
		supply[0] = +5.0
		supply[4] = -5.0

		self._add_edge_cost_capacity(G, capacities, 0, 1, 1.0, 3.0)
		self._add_edge_cost_capacity(G, capacities, 1, 4, 1.0, 3.0)
		self._add_edge_cost_capacity(G, capacities, 0, 2, 2.0, 3.0)
		self._add_edge_cost_capacity(G, capacities, 2, 4, 1.0, 3.0)

		solver = nk.flow.SuccessiveShortestPathMinCostFlow(G, self.CAPACITY, self.SUPPLY)
		solver.run()

		self.assertEqual(solver.getTotalCost(), 12.0)

		flow = solver.getFlow()
		self._check_flow_conservation(G, flow, supply)

	def testComplexMultiSourceNegativeCost(self):
		G = self._make_graph(6, weighted=True, directed=True, index_edges=True)
		capacities = G.attachEdgeAttribute(self.CAPACITY, float)
		supply = G.attachNodeAttribute(self.SUPPLY, float)

		for u in range(6):
			supply[u] = 0.0
		supply[0] = 4.0
		supply[1] = 3.0
		supply[5] = -7.0

		# edges with (cost, capacity)
		self._add_edge_cost_capacity(G, capacities, 0, 2,  2.0, 4.0)
		self._add_edge_cost_capacity(G, capacities, 2, 5,  1.0, 4.0)
		self._add_edge_cost_capacity(G, capacities, 0, 3,  1.0, 2.0)
		self._add_edge_cost_capacity(G, capacities, 3, 5,  3.0, 2.0)
		self._add_edge_cost_capacity(G, capacities, 1, 3, -1.0, 3.0)
		self._add_edge_cost_capacity(G, capacities, 3, 4,  1.0, 3.0)
		self._add_edge_cost_capacity(G, capacities, 4, 5,  2.0, 3.0)
		self._add_edge_cost_capacity(G, capacities, 1, 2,  3.0, 3.0)

		solver = nk.flow.SuccessiveShortestPathMinCostFlow(G, self.CAPACITY, self.SUPPLY)
		solver.run()

		self.assertAlmostEqual(solver.getTotalCost(), 18.0, places=12)

		flow = solver.getFlow()

		# Expected per-edge flows (others are 0.0)
		expected = {eid: 0.0 for eid in range(G.upperEdgeIdBound())}
		expected[G.edgeId(0, 2)] = 4.0
		expected[G.edgeId(2, 5)] = 4.0
		expected[G.edgeId(1, 3)] = 3.0
		expected[G.edgeId(3, 5)] = 2.0
		expected[G.edgeId(3, 4)] = 1.0
		expected[G.edgeId(4, 5)] = 1.0

		for (u, v) in G.iterEdges():
			eid = G.edgeId(u, v)
			self.assertAlmostEqual(flow[eid], expected[eid], places=12)

		self._check_flow_conservation(G, flow, supply)

	def testZeroCosts(self):
		G = self._make_graph(6, weighted=True, directed=True, index_edges=True)
		capacities = G.attachEdgeAttribute(self.CAPACITY, float)
		supply = G.attachNodeAttribute(self.SUPPLY, float)

		supply[0] = 4.0
		supply[1] = -2.0
		supply[2] = -2.0
		supply[3] = -4.0
		supply[4] = +2.0
		supply[5] = +2.0

		self._add_edge_cost_capacity(G, capacities, 0, 1, 0.0, 4.0)
		self._add_edge_cost_capacity(G, capacities, 0, 2, 0.0, 4.0)
		self._add_edge_cost_capacity(G, capacities, 4, 3, 0.0, 4.0)
		self._add_edge_cost_capacity(G, capacities, 5, 3, 0.0, 4.0)

		# capacity 0 edge (should carry no flow)
		self._add_edge_cost_capacity(G, capacities, 0, 3, 0.0, 0.0)

		solver = nk.flow.SuccessiveShortestPathMinCostFlow(G, self.CAPACITY, self.SUPPLY)
		solver.run()

		self.assertAlmostEqual(solver.getTotalCost(), 0.0, places=12)

		flow = solver.getFlow()

		expected = {eid: 0.0 for eid in range(G.upperEdgeIdBound())}
		expected[G.edgeId(0, 1)] = 2.0
		expected[G.edgeId(0, 2)] = 2.0
		expected[G.edgeId(4, 3)] = 2.0
		expected[G.edgeId(5, 3)] = 2.0
		# expected[G.edgeId(0, 3)] stays 0.0

		for (u, v) in G.iterEdges():
			eid = G.edgeId(u, v)
			self.assertAlmostEqual(flow[eid], expected[eid], places=12)

		self._check_flow_conservation(G, flow, supply)

if __name__ == "__main__":
	unittest.main()

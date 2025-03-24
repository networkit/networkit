#!/usr/bin/env python3
import unittest
import random
import networkit as nk
import pickle
import numpy as np
import scipy as sc

class TestGraph(unittest.TestCase):

	def testConstructorWithCooAndNumNodes(self):
		row = np.array([0, 1, 2])
		col = np.array([1, 2, 0])
		data = np.array([1.0, 2.0, 3.0])

		for weighted in [True, False]:
			for directed in [True, False]:
				G = nk.GraphFromCoo((data, (row, col)), n = 3, directed = directed, weighted = weighted)
				self.assertEqual(G.isDirected(), directed)
				self.assertEqual(G.isWeighted(), weighted)
				self.assertEqual(G.numberOfNodes(), 3)
				self.assertEqual(G.numberOfEdges(), 3)

	def testConstructorWithCoo(self):
		row = np.array([0, 1, 2])
		col = np.array([1, 2, 0])

		for weighted in [True, False]:
			for directed in [True, False]:
				G = nk.GraphFromCoo((row, col), directed = directed, weighted = weighted)
				self.assertEqual(G.isDirected(), directed)
				self.assertEqual(G.isWeighted(), weighted)
				self.assertEqual(G.numberOfNodes(), 3)
				self.assertEqual(G.numberOfEdges(), 3)

	def testConstructorWithMatrix(self):
		row = np.array([0, 1, 2])
		col = np.array([1, 2, 0])
		data = np.array([1.0, 2.0, 3.0])
		S = sc.sparse.coo_matrix((data, (row, col)), dtype = np.double)

		for weighted in [True, False]:
			for directed in [True, False]:
				G = nk.GraphFromCoo(S, n = 3, directed = directed, weighted = weighted)
				self.assertEqual(G.isDirected(), directed)
				self.assertEqual(G.isWeighted(), weighted)
				self.assertEqual(G.numberOfNodes(), 3)
				self.assertEqual(G.numberOfEdges(), 3)

	def testAddNodes(self):
		G = nk.Graph(0)
		G.addNodes(10)
		self.assertEqual(G.numberOfNodes(), 10)
		for u in range(G.numberOfNodes()):
			self.assertTrue(G.hasNode(u))

	def testAddEdgeWithMissingNodes(self):
		G = nk.Graph(0)

		with self.assertRaises(RuntimeError):
			G.addEdge(0, 1, addMissing = False)

		self.assertEqual(G.numberOfNodes(), 0)
		self.assertEqual(G.numberOfEdges(), 0)

		G.addEdge(0, 2, addMissing = True)

		self.assertTrue(G.hasNode(0))
		self.assertTrue(G.hasNode(2))
		self.assertTrue(G.hasEdge(0,2))
		self.assertFalse(G.hasNode(1)) #this node should not be active

		G.addEdge(0, 1, addMissing = True)

		self.assertEqual(G.numberOfNodes(), 3)
		self.assertEqual(G.numberOfEdges(), 2)

	def testEdgeIterator(self):
		for weighted in [True, False]:
			for directed in [True, False]:
				g = self.getSmallGraph(weighted, directed)
				for u, v in g.iterEdges():
					self.assertTrue(g.hasEdge(u, v))
				for u, v, w in g.iterEdgesWeights():
					self.assertTrue(g.hasEdge(u, v))
					self.assertEqual(g.weight(u, v), w)

	def testAddEdgesCooArrayWithArrayUndirected(self):
		# Using syntax similar to: https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_array.html
		row = np.array([0, 1, 2])
		col = np.array([1, 2, 0])
		data = np.array([1.0, 2.0, 3.0])

		G = nk.Graph(4)
		G.addEdges((data, (row, col)))

		self.assertEqual(G.numberOfEdges(), len(row))
		for i in range(len(row)):
			self.assertTrue(G.hasEdge(row[i], col[i]))
			self.assertTrue(G.hasEdge(col[i],row[i]))

	def testAddEdgesCooArrayWithArrayDirected(self):
		# Using syntax similar to: https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_array.html
		row = np.array([0, 1, 2])
		col = np.array([1, 2, 0])
		data = np.array([1.0, 2.0, 3.0])

		G = nk.Graph(4, directed = True)
		G.addEdges((data, (row, col)))

		self.assertEqual(G.numberOfEdges(), len(row))
		for i in range(len(row)):
			self.assertTrue(G.hasEdge(row[i], col[i]))
			self.assertFalse(G.hasEdge(col[i],row[i]))

	def testAddEdgesCooArrayWithArrayWeighted(self):
		# Using syntax similar to: https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_array.html
		row = np.array([0, 1, 2])
		col = np.array([1, 2, 0])
		data = np.array([1.0, 2.0, 3.0])

		G = nk.Graph(4, weighted = True)
		G.addEdges((data, (row, col)))

		for i in range(len(row)):
			self.assertEqual(G.weight(row[i], col[i]), data[i])

	def testAddEdgesCooArrayWithCoordsOnly(self):
		row = np.array([0, 1, 2])
		col = np.array([1, 2, 0])

		G = nk.Graph(4, weighted = True)
		G.addEdges((row, col))

		for i in range(len(row)):
			self.assertEqual(G.weight(row[i], col[i]), 1.0)

	def testAddEdgesCooArrayWithMatrix(self):

		S = sc.sparse.coo_matrix([[0, 1, 0],[1, 0, 0],[0.0, 0.0, 1.0]], dtype = np.double)

		G = nk.Graph(4)
		G.addEdges(S)

		self.assertEqual(G.numberOfEdges(), S.nnz)

		for i in range(len(S.row)):
			self.assertTrue(G.hasEdge(S.row[i], S.col[i]))

	def testGraphPickling(self):
		G = nk.Graph(2)
		G.addEdge(0,1)
		G.indexEdges()
		pickledGraph = pickle.dumps(G)
		G2 = pickle.loads(pickledGraph)
		self.assertEqual(G.numberOfNodes(), G2.numberOfNodes())
		self.assertEqual(G.numberOfEdges(), G2.numberOfEdges())
		self.assertEqual(G.edgeId(0,1), G2.edgeId(0,1))
	
	def testIterator(self):
		# Undirected
		G = nk.Graph(4, False, False)

		G.addEdge(0, 1)
		G.addEdge(0, 2)
		G.addEdge(1, 2)
		G.addEdge(3, 1)
		G.addEdge(3, 2)

		def nbrFunc(u, v, weight, edgeId):
			forEdgesNbrs.append(v)

		# Iterate through neighbours of node using iterNeighbors
		def nodeIter(node):
			nodeList = []
			nbrs = G.iterNeighbors(node)
			try:
				while nbrs is not None:
					nodeList.append(next(nbrs))
			except StopIteration:
				pass
			return nodeList

		for node in range(G.upperNodeIdBound()):
			forEdgesNbrs = []
			G.forEdgesOf(node, nbrFunc)
			nodeNbrs = nodeIter(node)
			self.assertEqual(sorted(forEdgesNbrs), sorted(nodeNbrs))

		# Directed
		G = nk.Graph(4, False, True)

		G.addEdge(0, 1)
		G.addEdge(0, 2)
		G.addEdge(3, 1)
		G.addEdge(3, 2)
		G.addEdge(1, 2)

		# Iterate through neighbours of node using iterNeighbors
		def nodeInIter(node):
			nodeList = []
			nbrs = G.iterInNeighbors(node)
			try:
				while nbrs is not None:
					nodeList.append(next(nbrs))
			except StopIteration:
				pass
			return nodeList

		for node in range(G.upperNodeIdBound()):
			forEdgesNbrs = []
			G.forInEdgesOf(node, nbrFunc)
			nodeInNbrs = nodeInIter(node)
			self.assertEqual(sorted(forEdgesNbrs), sorted(nodeInNbrs))				

	def testNeighbors(self):
		# Directed
		G = nk.Graph(4, False, True)

		G.addEdge(0, 1)
		G.addEdge(0, 2)
		G.addEdge(3, 1)
		G.addEdge(3, 2)
		G.addEdge(1, 2)

		def getNeighbors(u):
			neighbors = []
			for v in G.iterNeighbors(u):
				neighbors.append(v)
			return sorted(neighbors)

		def getInNeighbors(u):
			inNeighbors = []
			for v in G.iterInNeighbors(u):
				inNeighbors.append(v)
			return sorted(inNeighbors)

		self.assertListEqual(getNeighbors(0), [1, 2])
		self.assertListEqual(getNeighbors(1), [2])
		self.assertListEqual(getNeighbors(2), [])
		self.assertListEqual(getNeighbors(3), [1, 2])

		self.assertListEqual(getInNeighbors(0), [])
		self.assertListEqual(getInNeighbors(1), [0, 3])
		self.assertListEqual(getInNeighbors(2), [0, 1, 3])
		self.assertListEqual(getInNeighbors(3), [])

		# Undirected
		G = nk.Graph(4, False, False)

		G.addEdge(0, 1)
		G.addEdge(0, 2)
		G.addEdge(3, 1)
		G.addEdge(3, 2)
		G.addEdge(1, 2)

		self.assertEqual(getNeighbors(0), [1, 2])
		self.assertEqual(getNeighbors(1), [0, 2, 3])
		self.assertEqual(getNeighbors(2), [0, 1, 3])
		self.assertEqual(getNeighbors(3), [1, 2])

	def testNodeIterator(self):
		nk.setSeed(42, False)

		g = self.getSmallGraph()

		def doTest(g):
			nodes = []
			g.forNodes(lambda u: nodes.append(u))

			i = 0
			for u in g.iterNodes():
				self.assertEqual(u, nodes[i])
				i += 1

		doTest(g)
		g.removeNode(nk.graphtools.randomNode(g))
		g.removeNode(nk.graphtools.randomNode(g))
		doTest(g)		

	def testNodeAttributes(self):
		G = nk.Graph(5)

		for attType in [int, float, str]:
			attVals = [attType(i) for i in G.iterNodes()]

			attrs = G.attachNodeAttribute("attribute", attType)
			for u in G.iterNodes():
				attrs[u] = attVals[u]
				self.assertEqual(attrs[u], attVals[u])

			G.detachNodeAttribute("attribute")

	def testNodeAttributeReadWrite(self):
		G = nk.Graph(5)

		for attType in [int, float, str]:
			attVals = [attType(u) for u in G.iterNodes()]

			attrs = G.attachNodeAttribute("attribute", attType)
			attrsRead = G.attachNodeAttribute("attributeRead", attType)
			for u in G.iterNodes():
				attrs[u] = attVals[u]
			attrs.write('test.attrs')
			attrsRead.read('test.attrs')
			for u in G.iterNodes():
				self.assertEqual(attrs[u], attrsRead[u])

			G.detachNodeAttribute("attribute")
			G.detachNodeAttribute("attributeRead")

	def testEdgeAttributesMandatoryIndexing(self):
		G = nk.Graph(5)
		with self.assertRaises(Exception):
			G.attachEdgeAttribute("my Attr", int)

	def testEdgeAttributesByIndex(self):
		G = nk.Graph(5)
		G.indexEdges()

		G.addEdge(0, 1)
		G.addEdge(0, 2)
		G.addEdge(1, 3)
		G.addEdge(2, 4)

		for attType in [int, float, str]:
			attVals = [attType(u) for u in range(G.numberOfEdges())]

			attrs = G.attachEdgeAttribute("attribute", attType)
			for u in range(G.numberOfEdges()):
				attrs[u] = attVals[u]
				self.assertEqual(attrs[u], attVals[u])

			G.detachEdgeAttribute("attribute")

	def testEdgeAttributesByNodePair(self):
		G = nk.Graph(5)
		G.indexEdges()

		G.addEdge(0, 1)
		G.addEdge(0, 2)
		G.addEdge(1, 3)
		G.addEdge(2, 4)

		for attType in [int, float, str]:
			attVals = [attType(u+v) for u,v in G.iterEdges()]

			attrs = G.attachEdgeAttribute("attribute", attType)
			for u,v in G.iterEdges():
				attrs[u,v] = attVals[u]
				self.assertEqual(attrs[u,v], attVals[u])

			G.detachEdgeAttribute("attribute")
	
	def testEdgeAttributeReadWrite(self):
		G = nk.Graph(5)
		G.indexEdges()

		G.addEdge(0, 1)
		G.addEdge(0, 2)
		G.addEdge(1, 3)
		G.addEdge(2, 4)

		for attType in [int, float, str]:
			attVals = [attType(u+v) for u,v in G.iterEdges()]

			attrs = G.attachEdgeAttribute("attribute", attType)
			attrsRead = G.attachEdgeAttribute("attributeRead", attType)
			for u,v in G.iterEdges():
				attrs[u,v] = attVals[u]
			attrs.write('test.attrs')
			attrsRead.read('test.attrs')
			for u,v in G.iterEdges():
				self.assertEqual(attrs[u,v], attrsRead[u,v])

			G.detachEdgeAttribute("attribute")
			G.detachEdgeAttribute("attributeRead")

	def testEdgeAttributeGet(self):
		G = nk.Graph(5)
		G.indexEdges()

		G.addEdge(0, 1)
		G.addEdge(0, 2)
		G.addEdge(1, 3)
		G.addEdge(2, 4)

		for attType in [int, float, str]:
			attVals = [attType(u+v) for u,v in G.iterEdges()]

			attrs = G.attachEdgeAttribute("attribute", attType)
			for u,v in G.iterEdges():
				attrs[u,v] = attVals[u]
			
			attrsGet = G.getEdgeAttribute("attribute", attType)
			for u,v in G.iterEdges():
				self.assertEqual(attrs[u,v], attrsGet[u,v])

			G.detachEdgeAttribute("attribute")

	def testNodeAttributeReadWrite(self):
		G = nk.Graph(5)

		for attType in [int, float, str]:
			attVals = [attType(u) for u in G.iterNodes()]

			attrs = G.attachNodeAttribute("attribute", attType)
			for u in G.iterNodes():
				attrs[u] = attVals[u]

			attrsGet = G.getNodeAttribute("attribute", attType)
			for u in G.iterNodes():
				self.assertEqual(attrs[u], attrsGet[u])

			G.detachNodeAttribute("attribute")

	def testRandomEdgesReproducibility(self):
		numSamples = 10
		numSeeds = 3
		numRepeats = 10

		for directed in [False, True]:
			G = self.getSmallGraph(False, directed)

			results = [[] for i in range(numSeeds)]
			for repeats in range(numRepeats):
				for seed in range(numSeeds):
					nk.setSeed(seed, False)
					results[seed].append(nk.graphtools.randomEdges(G, numSamples))

			# assert results are different for different seeds
			for seed in range(1, numSeeds):
				self.assertNotEqual(results[0][0], results[seed][0])

			# assert results are consistent for same seeds
			for repeats in results:
				for x in repeats[1:]:
					self.assertEqual(repeats[0], x)

	def testRemoveAllEdges(self):
		for directed in [True, False]:
			for weighted in [True, False]:
				G = self.getSmallGraph(weighted, directed)
				G.removeAllEdges()
				self.assertEqual(G.numberOfEdges(), 0)
				G.forNodePairs(lambda u, v: self.assertFalse(G.hasEdge(u, v)))

	def testRemoveSelfLoops(self):
		for directed in [True, False]:
			for weighted in [True, False]:
				g  = self.getSmallGraph(weighted, directed)
				for i in range(10):
					u = nk.graphtools.randomNode(g)
					g.addEdge(u, u)

				nSelfLoops = g.numberOfSelfLoops()
				nEdges = g.numberOfEdges()

				g.removeSelfLoops()

				self.assertEqual(nEdges - nSelfLoops, g.numberOfEdges())
				self.assertEqual(g.numberOfSelfLoops(), 0)

				g.forNodes(lambda u: self.assertFalse(g.hasEdge(u, u)))

	def testRemoveMultiEdges(self):
		nMultiedges = 5
		for directed in [True, False]:
			for weighted in [True, False]:
				G = self.getSmallGraph(weighted, directed)
				nEdges = G.numberOfEdges()

				for _ in range(nMultiedges):
					u, v = nk.graphtools.randomEdge(G)
					G.addEdge(u, v)

				G.removeMultiEdges()
				self.assertEqual(G.numberOfEdges(), nEdges)

				for _ in range(nEdges):
					u, v = nk.graphtools.randomEdge(G)
					G.removeEdge(u, v)
					self.assertFalse(G.hasEdge(u, v))
				self.assertEqual(G.numberOfEdges(), 0)				

	def testSpanningForest(self):
		G = self.getSmallGraph()
		sf = nk.graph.SpanningForest(G)
		sf.run()
		F = sf.getForest()

		self.assertEqual(G.numberOfNodes(), F.numberOfNodes())
		for u in F.iterNodes():
			self.assertTrue(F.degree(u) > 0 or G.degree(u) == 0)

		for u, v in F.iterEdges():
			self.assertTrue(G.hasEdge(u, v))

	def getSmallGraph(self, weighted=False, directed=False):
		G = nk.Graph(4, weighted, directed)
		G.addEdge(0, 1, 1.0)
		G.addEdge(0, 2, 2.0)
		G.addEdge(3, 1, 4.0)
		G.addEdge(3, 2, 5.0)
		G.addEdge(1, 2, 3.0)

		return G
	
	def testEdgeSwap(self):
		G = nk.Graph(4, directed=True, weighted=True)
		G.addEdge(0,1, 1.0)
		G.addEdge(2,3, 2.0)
		G.indexEdges()
		G.swapEdge(0,1, 2,3) # now edges are:(0,3, 1,0) and (2,1, 2.0)
		self.assertTrue(G.hasEdge(0,3))
		self.assertTrue(G.hasEdge(2,1))
		# Test that the total edge weight remained the same.
		self.assertEqual(G.totalEdgeWeight(), 3.0)	

	def testIterNeighborsWeight(self):
		g = self.getSmallGraph(weighted=True, directed=True)
		for u in g.iterNodes():
			self.assertTrue(g.hasNode(u))
			for w in g.iterNeighborsWeights(u):
				self.assertEqual(w[1], g.weight(u, w[0]))
			for w in g.iterInNeighborsWeights(u):
				self.assertEqual(w[1], g.weight(w[0], u))		
	
	def testIncreaseEdgeWeight(self):
		G = nk.Graph(4, directed=True, weighted=True)
		G.addEdge(0,1, 1.0)
		G.increaseWeight(0,1, 2.0)
		self.assertEqual(G.weight(0,1), 3.0)
		G.sortEdges()
		G.compactEdges()
		self.assertTrue(G.checkConsistency())
	
	def testGraphEventOperators(self):
		up1 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.NODE_ADDITION, 0, 0, 0)	
		up2 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.NODE_ADDITION, 1, 1, 1)
		self.assertNotEqual(up1,up2)
		self.assertLess(up1,up2)
		self.assertLessEqual(up1,up2)
		self.assertGreater(up2,up1)
		self.assertGreaterEqual(up2,up1)	

	def testGraphFromGraphEventStream(self):
		up1 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.NODE_ADDITION, 0, 0, 0)	
		up2 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.NODE_ADDITION, 1, 1, 1)	
		up3 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.NODE_ADDITION, 2, 2, 2)		
		up4 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.NODE_ADDITION, 3, 3, 3)	
		upE1 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.EDGE_ADDITION, 1, 2, 2.0)	
		upE2 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.EDGE_ADDITION, 2, 3, 3.0)	
		upE3 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.EDGE_ADDITION, 1, 3, 1.0)	
		allUpdates=[up1,up2,up3,up4,upE1,upE2,upE3]
		G = nk.dynamics.graphFromStream(allUpdates, weighted=True,directed=True)
		self.assertEqual(G.numberOfNodes(), 4)
		self.assertEqual(G.numberOfEdges(), 3)		
	
	def testGraphDifference(self):
		G1 = self.getSmallGraph(weighted=True,directed=True)
		G2 = self.getSmallGraph(weighted=True,directed=True)
		update1 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.EDGE_ADDITION, 0, 3, 2.0)	
		update2 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.EDGE_ADDITION, 2, 3, 1.0)	
		updateG2 = nk.dynamics.GraphUpdater(G2)
		updateG2.update([update1,update2])
		GraphDiff = nk.dynamics.GraphDifference(G1,G2)
		GraphDiff.run()
		#testing graphDifference functions
		self.assertEqual(len(GraphDiff.getEdits()), 2)
		self.assertEqual(GraphDiff.getNumberOfEdits(), 2)
		self.assertEqual(GraphDiff.getNumberOfNodeAdditions(), 0)
		self.assertEqual(GraphDiff.getNumberOfNodeRestorations(), 0)
		self.assertEqual(GraphDiff.getNumberOfNodeRemovals(), 0)
		self.assertEqual(GraphDiff.getNumberOfEdgeAdditions(), 2)					
		self.assertEqual(GraphDiff.getNumberOfEdgeRemovals(), 0)
		self.assertEqual(GraphDiff.getNumberOfEdgeWeightUpdates(), 0)

	def testNodesInRandomOrder(self):
		G = nk.Graph(5)
		G.forNodesInRandomOrder(lambda u: self.assertTrue(G.hasNode(u)))

	def testNeighborsWeight(self):
		G = nk.Graph(4, directed=True, weighted=True)
		G.addEdge(0,2, 3.0)
		G.addEdge(1,0, 2.5)
		self.assertEqual(G.weightedDegree(0), 3.0)
		self.assertEqual(G.weightedDegreeIn(0), 2.5)
	
	def testRestoreNode(self):
		G=nk.Graph(2)
		G.removeNode(1)
		G.restoreNode(1)
		self.assertEqual(G.numberOfNodes(), 2)
		self.assertTrue(G.hasNode(1))	

	def testRandomMaximumSpanningForest(self):
		G = self.getSmallGraph()
		rmsf = nk.graph.RandomMaximumSpanningForest(G)
		rmsf.run()
		self.assertEqual(G.numberOfNodes(), rmsf.getMSF(True).numberOfNodes()) 

	def testUnionMaximumSpanningForest(self):
		G = self.getSmallGraph()
		un = nk.graph.UnionMaximumSpanningForest(G)
		un.run()
		self.assertEqual(un.getUMSF().numberOfNodes(),4)
		self.assertEqual(un.getUMSF().numberOfEdges(),5)
		self.assertTrue(un.inUMST(0,1))	

if __name__ == "__main__":
	unittest.main()
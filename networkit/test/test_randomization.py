#!/usr/bin/env python3
import unittest
import sys
from copy import copy
import networkit as nk

def check_graphs(G1, G2):
    assert(G1.numberOfNodes() == G2.numberOfNodes())
    assert(G1.numberOfEdges() == G2.numberOfEdges())

    failed = False
    for i in range(G1.numberOfNodes()):
        if (G1.degree(i) != G2.degree(i)):
            print("Degree mismatch of node %d (%d != %d)" % (i, G1.degree(i), G2.degree(i)))
            failed = True

        if (G1.degreeIn(i) != G2.degreeIn(i)):
            print("In-Degree mismatch of node %d (%d != %d)" % (i, G1.degree(i), G2.degree(i)))
            failed = True

    if failed:
        nk.overview(G1)
        nk.overview(G2)
        raise RuntimeError("Degree mismatch")

class TestRandomization(unittest.TestCase):
    def setUp(self):
        self.graphs = []
        self.graphs.append(nk.generators.HyperbolicGenerator(1001,  50).generate())
        self.graphs.append(nk.generators.HyperbolicGenerator(1002,   5).generate())
        self.graphs.append(nk.generators.HyperbolicGenerator(10003, 10).generate())
        self.graphs.append(nk.generators.ErdosRenyiGenerator(1004, 0.005).generate())
        self.graphs.append(nk.generators.ErdosRenyiGenerator(1005, 0.05).generate())
        self.graphs.append(nk.generators.ErdosRenyiGenerator(1004, 0.005, True).generate())
        self.graphs.append(nk.generators.ErdosRenyiGenerator(1005, 0.05, True).generate())

    def testGlobalCurveball(self):
        for G in self.graphs:
            algo = nk.randomization.GlobalCurveball(G, 5)
            algo.run()
            G2 = algo.getGraph()
            check_graphs(G, G2)

    def testGlobalCurveballWithSelfloops(self):
        for G in self.graphs:
            if not G.isDirected(): continue

            algo = nk.randomization.GlobalCurveball(G, 5, True, True)
            algo.run()
            G2 = algo.getGraph()
            check_graphs(G, G2)

    def testGlobalCurveballWithPreprocessing(self):
        for G in self.graphs:
            algo = nk.randomization.GlobalCurveball(G, 5, False, True)
            algo.run()
            G2 = algo.getGraph()
            check_graphs(G, G2)

    def testCurveballWithGlobal(self):
        for G in self.graphs:
            if G.isDirected(): continue

            n = G.numberOfNodes()
            ts = nk.randomization.CurveballGlobalTradeGenerator(5, n).generate()
            algo = nk.randomization.Curveball(G)
            algo.run(ts)
            algo.run(ts)
            G2 = algo.getGraph()
            check_graphs(G, G2)

    def testCurveballWithUniform(self):
        for G in self.graphs:
            if G.isDirected(): continue

            n = G.numberOfNodes()
            ts = nk.randomization.CurveballUniformTradeGenerator(5 * n, n).generate()
            algo = nk.randomization.Curveball(G)
            algo.run(ts)
            algo.run(ts)
            G2 = algo.getGraph()
            check_graphs(G, G2)

    def testDegreePreservingShuffle(self):
        for G in self.graphs:
            dps = nk.randomization.DegreePreservingShuffle(G)
            dps.run()
            G2 = dps.getGraph()
            check_graphs(G, G2)
            perm = dps.getPermutation()
            for u in G.iterNodes():
                self.assertEqual(G.degree(u), G.degree(perm[u]))
                self.assertEqual(G.degreeIn(u), G.degreeIn(perm[u]))

    def testDegreePreservingShuffleDirectedTriangle(self):
        """Test whether a directed triangle is reoriented in 50% of cases"""
        G = nk.Graph(3, False, True)
        G.addEdge(0, 1)
        G.addEdge(1, 2)
        G.addEdge(2, 0)

        num_clockwise = 0
        num_iterations = 1000

        for i in range(num_iterations):
            dps = nk.randomization.DegreePreservingShuffle(G)
            dps.run()
            G2 = dps.getGraph()
            check_graphs(G, G2)

            # check orientation
            clockwise = G2.hasEdge(0, 1)
            anticlkw  = G2.hasEdge(1, 0)
            self.assertNotEqual(clockwise, anticlkw)

            num_clockwise += clockwise
            G = G2

        # confidence interval with an error rate of ~ 1e-6
        self.assertGreater(num_clockwise, 400)
        self.assertLess   (num_clockwise, 600)

    def testEdgeSwitching(self):
        for numSwitchesPerEdge, preShuffle in [(0, False), (0, True), (1, False)]:
            G = nk.generators.ErdosRenyiGenerator(100, 0.1).generate()
            algo = nk.randomization.EdgeSwitching(G, numSwitchesPerEdge + 1, preShuffle)
            self.assertEqual(algo.getNumberOfSwitchesPerEdge(), numSwitchesPerEdge + 1)
            algo.setNumberOfSwitchesPerEdge(numSwitchesPerEdge)
            self.assertEqual(algo.getNumberOfSwitchesPerEdge(), numSwitchesPerEdge)
            numSwitches = numSwitchesPerEdge * G.numberOfEdges()

            algo.run()
            G1 = algo.getGraph()

            self.assertEqual(G1.numberOfNodes(), G.numberOfNodes())
            self.assertEqual(G1.numberOfEdges(), G.numberOfEdges())
            if numSwitches > 0 or preShuffle:
                self.assertNotEqual(sorted(G.iterEdges()), sorted(G1.iterEdges()))
            else:
                self.assertEqual(sorted(G.iterEdges()), sorted(G1.iterEdges()))

            self.assertGreaterEqual(algo.getNumberOfAffectedEdges(), numSwitches // 2)

    def testEdgeSwitchingInplace(self):
        for numSwitchesPerEdge in [0, 1]:
            G = nk.generators.ErdosRenyiGenerator(100, 0.1).generate()
            G_old = copy(G)
            algo = nk.randomization.EdgeSwitchingInPlace(G, numSwitchesPerEdge + 1)
            self.assertEqual(algo.getNumberOfSwitchesPerEdge(), numSwitchesPerEdge + 1)
            algo.setNumberOfSwitchesPerEdge(numSwitchesPerEdge)
            self.assertEqual(algo.getNumberOfSwitchesPerEdge(), numSwitchesPerEdge)
            numSwitches = numSwitchesPerEdge * G.numberOfEdges()

            algo.run()

            self.assertEqual(G_old.numberOfNodes(), G.numberOfNodes())
            self.assertEqual(G_old.numberOfEdges(), G.numberOfEdges())

            if numSwitches > 0:
                self.assertNotEqual(sorted(G_old.iterEdges()), sorted(G.iterEdges()))
            else:
                self.assertEqual(sorted(G_old.iterEdges()), sorted(G.iterEdges()))

            self.assertGreaterEqual(algo.getNumberOfAffectedEdges(), numSwitches // 2)

    def testEdgeSwitchingInplaceRefCount(self):
        G = nk.generators.ErdosRenyiGenerator(10, 0.1).generate()
        rc_initial = sys.getrefcount(G)
        algo = nk.randomization.EdgeSwitchingInPlace(G)
        self.assertGreater(sys.getrefcount(G), rc_initial)
        algo.run()
        self.assertGreater(sys.getrefcount(G), rc_initial)
        del(algo)
        self.assertEqual(sys.getrefcount(G), rc_initial)


if __name__ == "__main__":
    unittest.main()

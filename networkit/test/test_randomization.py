#!/usr/bin/env python3
import unittest
import networkit as nk

def check_graphs(G1, G2):
    assert(G1.numberOfNodes() == G2.numberOfNodes())
    assert(G1.numberOfEdges() == G2.numberOfEdges())

    failed = False
    for i in range(G1.numberOfNodes()):
        if (G1.degree(i) != G2.degree(i)):
            print("Degree mismatch of node %d (%d != %d)" % (i, G1.degree(i), G2.degree(i)))
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

    def test_global_curveball(self):
        for G in self.graphs:
            algo = nk.randomization.GlobalCurveball(G, 5)
            algo.run()
            G2 = algo.getGraph()
            check_graphs(G, G2)

    def test_curveball_with_global(self):
        for G in self.graphs:
            n = G.numberOfNodes()
            ts = nk.randomization.CurveballGlobalTradeGenerator(5, n).generate()
            algo = nk.randomization.Curveball(G)
            algo.run(ts)
            algo.run(ts)
            G2 = algo.getGraph()
            check_graphs(G, G2)

    def test_curveball_with_uniform(self):
        for G in self.graphs:
            n = G.numberOfNodes()
            ts = nk.randomization.CurveballUniformTradeGenerator(5 * n, n).generate()
            algo = nk.randomization.Curveball(G)
            algo.run(ts)
            algo.run(ts)
            G2 = algo.getGraph()
            check_graphs(G, G2)

if __name__ == "__main__":
    unittest.main()

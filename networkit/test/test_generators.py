#!/usr/bin/env python3
import unittest
import networkit as nk

class TestGraph(unittest.TestCase):
    def testPubWebGenerator(self):
        n = 100
        gen = nk.generators.PubWebGenerator(n, 8, 0.1, 4)
        G = gen.generate()
        coords = gen.getCoordinates()

        self.assertEqual(G.numberOfNodes(), n)
        self.assertEqual(len(coords), n)
        self.assertTrue(all(len(pt) == 2 for pt in coords))

        # check that all points are within the unit square and that they are non-empt
        self.assertTrue(all(0.0 <= pt[0] <= 1 for pt in coords))
        self.assertTrue(any(0.0 <  pt[0] <  1 for pt in coords))
        self.assertTrue(all(0.0 <= pt[1] <= 1 for pt in coords))
        self.assertTrue(any(0.0 <  pt[1] <  1 for pt in coords))

    def testDynamicPubWebGenerator(self):
        n = 100
        gen = nk.generators.DynamicPubWebGenerator(n, 8, 0.1, 4)
        G = gen.getGraph()
        coords = gen.getCoordinates()

        self.assertEqual(G.numberOfNodes(), n)
        self.assertEqual(len(coords), n)
        self.assertTrue(all(len(pt) == 2 for pt in coords))

        # check that all points are within the unit square and that they are non-empt
        self.assertTrue(all(0.0 <= pt[0] <= 1 for pt in coords))
        self.assertTrue(any(0.0 <  pt[0] <  1 for pt in coords))
        self.assertTrue(all(0.0 <= pt[1] <= 1 for pt in coords))
        self.assertTrue(any(0.0 <  pt[1] <  1 for pt in coords))

        gen.generate(10)

        newCoords = gen.getNewCoordinates()
        self.assertGreater(len(newCoords), 0)
        self.assertEqual(newCoords[0][0], n)

    def testWattsStrogatzGenerator(self):
        n = 12
        nbrs = 5
        G = nk.generators.WattsStrogatzGenerator(n, nbrs, 0.5).generate()

        self.assertEqual(G.numberOfNodes(), n)
        self.assertEqual(G.numberOfEdges(), n*nbrs)

    def testClusteredRandomGraphGenerator(self):
        n, k = 100, 10
        pIntra, pInter = 0.1, 0.05
        gen = nk.generators.ClusteredRandomGraphGenerator(n, k, pIntra, pInter)
        G = gen.generate()
        self.assertEqual(G.numberOfNodes(), n)
        nCommunities = len(set(gen.getCommunities().getVector()))
        self.assertGreater(nCommunities, 0)
        self.assertLessEqual(nCommunities, k)


    def testEdgeSwitchingMarkovChainGenerator(self):
        n = 100
        while True:
            degseq = nk.generators.PowerlawDegreeSequence(2, n).run().getDegreeSequence(n)
            gen = nk.generators.EdgeSwitchingMarkovChainGenerator(degseq, False, 10)
            if not gen.isRealizable():
                continue

            G = gen.generate()
            self.assertEqual(G.numberOfNodes(), n)
            for u in range(G.numberOfNodes()):
                self.assertEqual(G.degree(u), degseq[u])

            break

if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
import unittest
import networkit as nk


class TestGraph(unittest.TestCase):

    def setUp(self):
        self.G = nk.readGraph("input/looptest1.gml", nk.Format.GML)

    def testClusteredRandomGraphGenerator(self):
        n, k = 100, 10
        pIntra, pInter = 0.1, 0.05
        gen = nk.generators.ClusteredRandomGraphGenerator(n, k, pIntra, pInter)
        G = gen.generate()
        self.assertEqual(G.numberOfNodes(), n)
        nCommunities = len(set(gen.getCommunities().getVector()))
        self.assertGreater(nCommunities, 0)
        self.assertLessEqual(nCommunities, k)

    def testBarabasiAlbertGenerator(self):
        g = nk.generators.BarabasiAlbertGenerator(2, 50, 1)
        g.fit(self.G)
        self.assertEqual(g.generate().numberOfNodes(), 50)
        self.assertEqual(g.generate().numberOfEdges(), 98)

    def testBTERReplicator(self):
        g = nk.generators.BTERReplicator(self.G)
        g.setPaths("./BTER.data")
        g.fit(self.G)

    def testChungLuGenerator(self):
        g = nk.generators.ChungLuGenerator([1.3, 1.5, 1.7])
        g.fit(self.G)
        self.assertEqual(g.generate().numberOfNodes(), 3)

    def testDorogovtsevMendesGenerator(self):
        g = nk.generators.DorogovtsevMendesGenerator(100)
        g.fit(self.G)
        self.assertEqual(g.generate().numberOfNodes(), 100)
        self.assertEqual(g.generate().numberOfEdges(), 197)

    def testDynamicHyperbolicGenerator(self):
        g = nk.generators.DynamicHyperbolicGenerator(100)
        g.generate(5)
        res = g.getGraph()
        coords = g.getCoordinates()
        self.assertEqual(res.numberOfNodes(), 100)
        self.assertEqual(len(coords), 100)

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
        self.assertTrue(any(0.0 < pt[0] < 1 for pt in coords))
        self.assertTrue(all(0.0 <= pt[1] <= 1 for pt in coords))
        self.assertTrue(any(0.0 < pt[1] < 1 for pt in coords))

        gen.generate(10)

        newCoords = gen.getNewCoordinates()
        self.assertGreater(len(newCoords), 0)
        self.assertEqual(newCoords[0][0], n)

    def testDynamicForestFireGenerator(self):
        dyn = nk.generators.DynamicForestFireGenerator(0.5, directed=True)
        G = dyn.generate(10)
        self.assertGreaterEqual(len(G), 20)

    def testEdgeSwitchingMarkovChainGenerator(self):
        n = 100
        while True:
            degseq = (
                nk.generators.PowerlawDegreeSequence(2, n).run().getDegreeSequence(n)
            )
            gen = nk.generators.EdgeSwitchingMarkovChainGenerator(degseq, False, 10)
            if not gen.isRealizable():
                continue

            G = gen.generate()
            self.assertEqual(G.numberOfNodes(), n)
            for u in range(G.numberOfNodes()):
                self.assertEqual(G.degree(u), degseq[u])

            break

    def testErdosRenyiGenerator(self):
        g = nk.generators.ErdosRenyiGenerator(10, 0.15, True)
        g.fit(self.G)
        self.assertEqual(g.generate().numberOfNodes(), 10)

    def testHavelHakimiGenerator(self):
        g = nk.generators.HavelHakimiGenerator([2, 2, 2, 2])
        self.assertTrue(g.isRealizable())
        self.assertTrue(g.getRealizable())

    def testHavelHakimiGeneratorFit(self):
        g = nk.generators.HavelHakimiGenerator([2, 2, 2, 2])
        g.fit(self.G)
        self.assertEqual(g.generate().numberOfNodes(), 4)
        self.assertEqual(g.generate().numberOfEdges(), 4)

    def testHyperbolicGenerator(self):
        g = nk.generators.HyperbolicGenerator(100)
        g.setLeafCapacity(5)
        g.setBalance(0.5)
        g.setTheoreticalSplit(False)
        self.assertEqual(g.generate().numberOfNodes(), 100)

    def testLFRGenerator(self):
        g = nk.generators.LFRGenerator(5)
        g.setDegreeSequence([0, 0, 0, 0, 0])
        g.setMu(1.5)
        g.setCommunitySizeSequence([1, 1, 1, 1, 1])
        self.assertEqual(g.generate().numberOfNodes(), 5)
        g.run()
        self.assertEqual(g.getGraph().numberOfEdges(), 0)
        self.assertEqual(g.getPartition().numberOfElements(), 5)

    def testLFRGeneratorFit(self):
        g = nk.generators.LFRGenerator(5)
        g.setDegreeSequence([0, 0, 0, 0, 0])
        g.setMu(1.5)
        g.setCommunitySizeSequence([1, 1, 1, 1, 1])
        g.fit(self.G)
        self.assertEqual(g.generate().numberOfNodes(), 5)
        self.assertEqual(g.generate().numberOfEdges(), 0)

    def testLFRGeneratorPowerlawDegreeSequence(self):
        # Example with generatingPowerlaw functions
        g2 = nk.generators.LFRGenerator(100)
        g2.generatePowerlawDegreeSequence(10, 20, -2)
        g2.generatePowerlawCommunitySizeSequence(10, 50, -1)
        g2.setMu(0.5)
        self.assertEqual(g2.generate().numberOfNodes(), 100)
        g2.run()
        self.assertEqual(g2.getGraph().numberOfNodes(), 100)
        self.assertEqual(g2.getPartition().numberOfElements(), 100)

    def testLFRGeneratorPowerlawDegreeSequenceFit(self):
        g2 = nk.generators.LFRGenerator(10)
        g2.generatePowerlawDegreeSequence(5, 6, -2)
        g2.generatePowerlawCommunitySizeSequence(5, 6, -1)
        g2.setMu(0.5)
        g2.fit(self.G, vanilla=False, scale=2)
        self.assertEqual(g2.generate().numberOfNodes(), 10)

    def testLFRGeneratorPowerlawDegreeSequenceFitVanilla(self):
        g2 = nk.generators.LFRGenerator(10)
        g2.generatePowerlawDegreeSequence(5, 6, -2)
        g2.generatePowerlawCommunitySizeSequence(5, 6, -1)
        g2.setMu(0.5)
        g2.fit(self.G, vanilla=True, plfit=False)
        self.assertEqual(g2.generate().numberOfNodes(), 10)

    def testPowerLawDegreeSequence(self):
        n = 100
        PowGen = nk.generators.PowerlawDegreeSequence(2, n)
        PowGen.setMinimumFromAverageDegree(12.0)
        PowGen.setGammaFromAverageDegree(3.0)
        PowGen.setGamma(5.0)
        PowGen.run()
        self.assertAlmostEqual(PowGen.getExpectedAverageDegree(), 86.139, 3)
        self.assertEqual(PowGen.getMinimumDegree(), 4)
        self.assertEqual(PowGen.getMaximumDegree(), 100)
        self.assertEqual(PowGen.getGamma(), 5.0)
        # Algorithm is not deterministic
        self.assertIsInstance(PowGen.getDegree(), int)

    def testPowerLawDegreeSequenceFromGraph(self):
        PowGen = nk.generators.PowerlawDegreeSequence(self.G)
        PowGen.run()

        degG = nk.centrality.DegreeCentrality(self.G).run().scores()
        minDegG = min(degG)
        maxDegG = max(degG)
        self.assertEqual(PowGen.getMinimumDegree(), minDegG)
        self.assertEqual(PowGen.getMaximumDegree(), maxDegG)

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
        self.assertTrue(any(0.0 < pt[0] < 1 for pt in coords))
        self.assertTrue(all(0.0 <= pt[1] <= 1 for pt in coords))
        self.assertTrue(any(0.0 < pt[1] < 1 for pt in coords))

    def testRmatGenerator(self):
        g = nk.generators.RmatGenerator(2, 1, 0.25, 0.25, 0.25, 0.25)
        g.setPaths("./rMat.data")
        resGraph = g.generate()
        self.assertAlmostEqual(resGraph.numberOfNodes(), 4)
        self.assertAlmostEqual(resGraph.numberOfEdges(), 4)

    def testRmatGeneratorFit(self):
        g = nk.generators.RmatGenerator(2, 1, 0.25, 0.25, 0.25, 0.25)
        g.setPaths("./rMat.data")
        inputGraph = nk.Graph(3)
        inputGraph.addEdge(0, 1)
        inputGraph.addEdge(1, 2)
        g.fit(inputGraph, kronfit=False)
        self.assertEqual(g.generate().numberOfNodes(), 4)
        self.assertEqual(g.generate().numberOfEdges(), 4)

    def testWattsStrogatzGenerator(self):
        n = 12
        nbrs = 5
        G = nk.generators.WattsStrogatzGenerator(n, nbrs, 0.5).generate()

        self.assertEqual(G.numberOfNodes(), n)
        self.assertEqual(G.numberOfEdges(), n * nbrs)


if __name__ == "__main__":
    unittest.main()

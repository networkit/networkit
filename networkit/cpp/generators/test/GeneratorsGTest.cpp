/*
 *  GeneratorsTest.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <numeric>

#include <networkit/generators/BarabasiAlbertGenerator.hpp>
#include <networkit/generators/ChungLuGenerator.hpp>
#include <networkit/generators/ChungLuGeneratorAlamEtAl.hpp>
#include <networkit/generators/ClusteredRandomGraphGenerator.hpp>
#include <networkit/generators/ConfigurationModel.hpp>
#include <networkit/generators/DorogovtsevMendesGenerator.hpp>
#include <networkit/generators/DynamicBarabasiAlbertGenerator.hpp>
#include <networkit/generators/DynamicDorogovtsevMendesGenerator.hpp>
#include <networkit/generators/DynamicForestFireGenerator.hpp>
#include <networkit/generators/DynamicGraphSource.hpp>
#include <networkit/generators/DynamicHyperbolicGenerator.hpp>
#include <networkit/generators/DynamicPathGenerator.hpp>
#include <networkit/generators/DynamicPubWebGenerator.hpp>
#include <networkit/generators/EdgeSwitchingMarkovChainGenerator.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/generators/HavelHakimiGenerator.hpp>
#include <networkit/generators/HyperbolicGenerator.hpp>
#include <networkit/generators/LFRGenerator.hpp>
#include <networkit/generators/MocnikGenerator.hpp>
#include <networkit/generators/MocnikGeneratorBasic.hpp>
#include <networkit/generators/PowerlawDegreeSequence.hpp>
#include <networkit/generators/PubWebGenerator.hpp>
#include <networkit/generators/RegularRingLatticeGenerator.hpp>
#include <networkit/generators/RmatGenerator.hpp>
#include <networkit/generators/StochasticBlockmodel.hpp>
#include <networkit/generators/WattsStrogatzGenerator.hpp>

#include <networkit/auxiliary/MissingMath.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/community/ClusteringGenerator.hpp>
#include <networkit/community/Modularity.hpp>
#include <networkit/community/PLM.hpp>
#include <networkit/community/PLP.hpp>
#include <networkit/dynamics/GraphUpdater.hpp>
#include <networkit/global/ClusteringCoefficient.hpp>
#include <networkit/io/DotGraphWriter.hpp>
#include <networkit/io/GraphIO.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/METISGraphWriter.hpp>
#include <networkit/viz/PostscriptWriter.hpp>

namespace NetworKit {

class GeneratorsGTest : public testing::Test {

protected:
    template <class Type>
    std::vector<Type> readVector(std::string_view path) const {
        std::ifstream inputFile(path.data());
        Type cur;
        std::vector<Type> data;

        while (inputFile >> cur)
            data.push_back(cur);

        return data;
    }

public:
    vector<double> getAngles(const DynamicHyperbolicGenerator &dynGen) { return dynGen.angles; }

    vector<double> getRadii(const DynamicHyperbolicGenerator &dynGen) { return dynGen.radii; }
};

TEST_F(GeneratorsGTest, testClusteredRandomGraphGenerator) {
    Aux::Random::setSeed(42, false);
    const count n = 100, c = 10;
    const double pin = 0.5, pout = 0.01;
    ClusteredRandomGraphGenerator gen(n, c, pin, pout);
    Graph G = gen.generate();
    Partition part = gen.getCommunities();
    count nCommunities = part.getSubsetIds().size();
    EXPECT_EQ(n, G.numberOfNodes());
    EXPECT_EQ(n, G.upperNodeIdBound());
    EXPECT_GE(nCommunities, 1);
    EXPECT_LE(nCommunities, c);
}

TEST_F(GeneratorsGTest, testClusteredRandomGraphGeneratorCompleteIsolatedCommunities) {
    const count n = 100, c = 10;
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        ClusteredRandomGraphGenerator gen(n, c, 1, 0);
        auto G = gen.generate();
        EXPECT_EQ(n, G.numberOfNodes());
        EXPECT_EQ(n, G.upperNodeIdBound());
        EXPECT_EQ(G.numberOfSelfLoops(), 0);
        { // Check no multiple edges
            auto G1 = G;
            G1.removeMultiEdges();
            EXPECT_EQ(G.numberOfEdges(), G1.numberOfEdges());
        }

        const auto commPartition = gen.getCommunities();
        const auto communities = commPartition.getSubsets();
        EXPECT_GE(communities.size(), 1);
        EXPECT_LE(communities.size(), c);

        // Check that each community is an isolated clique
        const auto subsetSizes = commPartition.subsetSizeMap();
        for (const auto &community : communities) {
            for (node u : community) {
                const index uIdx = commPartition.subsetOf(u);
                EXPECT_EQ(G.degree(u), subsetSizes.at(uIdx) - 1);
                G.forNeighborsOf(u, [&](node v) { EXPECT_EQ(uIdx, commPartition.subsetOf(v)); });
            }
        }
    }
}

TEST_F(GeneratorsGTest, testClusteredRandomGraphGeneratorCompleteCommunities) {
    const count n = 100, c = 10;
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        ClusteredRandomGraphGenerator gen(n, c, 1, 0.01);
        auto G = gen.generate();
        EXPECT_EQ(n, G.numberOfNodes());
        EXPECT_EQ(n, G.upperNodeIdBound());
        EXPECT_EQ(G.numberOfSelfLoops(), 0);
        { // Check no multiple edges
            auto G1 = G;
            G1.removeMultiEdges();
            EXPECT_EQ(G.numberOfEdges(), G1.numberOfEdges());
        }

        const auto commPartition = gen.getCommunities();
        const auto communities = commPartition.getSubsets();
        EXPECT_GE(communities.size(), 1);
        EXPECT_LE(communities.size(), c);

        // Check that each community is a clique
        const auto subsetSizes = commPartition.subsetSizeMap();
        for (const auto &community : communities) {
            for (node u : community) {
                const index uIdx = commPartition.subsetOf(u);
                EXPECT_GE(G.degree(u), subsetSizes.at(uIdx) - 1);
                count leftToVisit = community.size() - 1;
                G.forNeighborsOf(u, [&](node v) {
                    if (uIdx == commPartition.subsetOf(v))
                        --leftToVisit;
                });
                EXPECT_EQ(leftToVisit, 0);
            }
        }
    }
}

TEST_F(GeneratorsGTest, testClusteredRandomGraphGeneratorIndependentCommunities) {
    const count n = 100, c = 10;
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        ClusteredRandomGraphGenerator gen(n, c, 0, 1);
        auto G = gen.generate();
        EXPECT_EQ(n, G.numberOfNodes());
        EXPECT_EQ(n, G.upperNodeIdBound());
        EXPECT_EQ(G.numberOfSelfLoops(), 0);
        { // Check no multiple edges
            auto G1 = G;
            G1.removeMultiEdges();
            EXPECT_EQ(G.numberOfEdges(), G1.numberOfEdges());
        }

        const auto commPartition = gen.getCommunities();
        const auto communities = commPartition.getSubsets();
        EXPECT_GE(communities.size(), 1);
        EXPECT_LE(communities.size(), c);

        // Check that each community is an independent set, and the vertices in
        // each community are connected to all the vertices outside the community.
        const auto subsetSizes = commPartition.subsetSizeMap();
        for (const auto &community : communities) {
            for (node u : community) {
                EXPECT_EQ(G.degree(u), n - community.size());
                G.forNeighborsOf(u, [&](node v) {
                    EXPECT_NE(commPartition.subsetOf(u), commPartition.subsetOf(v));
                });
            }
        }
    }
}

TEST_F(GeneratorsGTest, testClusteredRandomGraphGeneratorSparseCommunities) {
    const count n = 100, c = 10;
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        ClusteredRandomGraphGenerator gen(n, c, 0.01, 1);
        auto G = gen.generate();
        EXPECT_EQ(n, G.numberOfNodes());
        EXPECT_EQ(n, G.upperNodeIdBound());
        EXPECT_EQ(G.numberOfSelfLoops(), 0);
        { // Check no multiple edges
            auto G1 = G;
            G1.removeMultiEdges();
            EXPECT_EQ(G.numberOfEdges(), G1.numberOfEdges());
        }

        const auto commPartition = gen.getCommunities();
        const auto communities = commPartition.getSubsets();
        EXPECT_GE(communities.size(), 1);
        EXPECT_LE(communities.size(), c);

        // Check that the vertices in each community are connected to all the
        // vertices outside the community.
        const auto subsetSizes = commPartition.subsetSizeMap();
        for (const auto &community : communities) {
            for (node u : community) {
                EXPECT_GE(G.degree(u), n - community.size());
                count leftToVisit = n - community.size();
                G.forNeighborsOf(u, [&](node v) {
                    if (commPartition.subsetOf(u) != commPartition.subsetOf(v))
                        --leftToVisit;
                });
                EXPECT_EQ(leftToVisit, 0);
            }
        }
    }
}

TEST_F(GeneratorsGTest, testConfigurationModelGeneratorWithEdgeSwitchingOnRealSequence) {
    METISGraphReader reader;
    auto graphs = {"input/jazz.graph", "input/lesmis.graph"};

    for (auto path : graphs) {
        Graph G = reader.read(path);
        count n = G.upperNodeIdBound();
        std::vector<count> sequence(n);
        G.forNodes([&](node u) { sequence[u] = G.degree(u); });

        bool skipTest = false;
        EdgeSwitchingMarkovChainGenerator gen(sequence, skipTest);
        Graph G2 = gen.generate();

        count volume = std::accumulate(sequence.begin(), sequence.end(), 0);
        EXPECT_EQ(volume, 2 * G2.numberOfEdges());

        std::vector<count> testSequence(n);
        G2.forNodes([&](node u) { testSequence[u] = G2.degree(u); });
        Aux::Parallel::sort(testSequence.begin(), testSequence.end(), std::greater<count>());
        Aux::Parallel::sort(sequence.begin(), sequence.end(), std::greater<count>());

        for (index i = 0; i < n; ++i) {
            EXPECT_EQ(sequence[i], testSequence[i]);
        }
    }
}

TEST_F(GeneratorsGTest, testConfigurationModelGeneratorWithRejectionSamplingOnRealSequence) {
    METISGraphReader reader;
    Aux::Random::setSeed(42, false);
    auto graphs = {"input/jazz.graph", "input/lesmis.graph"};

    for (auto path : graphs) {
        Graph G = reader.read(path);
        count n = G.numberOfNodes();
        std::vector<count> sequence(n);
        G.forNodes([&](node u) { sequence[u] = G.degree(u); });

        ConfigurationModel gen(sequence);
        Graph G2 = gen.generate();

        count volume = std::accumulate(sequence.begin(), sequence.end(), 0);
        EXPECT_EQ(volume, 2 * G2.numberOfEdges());

        std::vector<count> testSequence(n);
        G2.forNodes([&](node u) { testSequence[u] = G2.degree(u); });
        Aux::Parallel::sort(testSequence.begin(), testSequence.end(), std::greater<count>());
        Aux::Parallel::sort(sequence.begin(), sequence.end(), std::greater<count>());

        for (index i = 0; i < n; ++i) {
            EXPECT_EQ(sequence[i], testSequence[i]);
        }
    }
}

TEST_F(GeneratorsGTest, testDynamicBarabasiAlbertGeneratorSingleStep) {
    count k = 2; // number of edges added per node
    DynamicGraphSource *gen = new DynamicBarabasiAlbertGenerator(k);
    GraphEventProxy *Gproxy = gen->newGraph();
    Graph *G = Gproxy->G;

    gen->initializeGraph();

    count nPre = G->numberOfNodes();
    count mPre = G->numberOfEdges();
    EXPECT_EQ(k, nPre) << "graph should have been initialized to k nodes";
    EXPECT_EQ(k - 1, mPre)
        << "graph should have been initialized to a path of k nodes which means k-1 edges";

    // perform single preferential attachment step
    gen->generate();

    count nPost = G->numberOfNodes();
    count mPost = G->numberOfEdges();
    EXPECT_EQ(nPre + 1, nPost) << "one more node should have been added";
    EXPECT_EQ(mPre + k, mPost) << "k edges should have been added";

    delete gen;
    delete Gproxy;
    delete G;
}

TEST_F(GeneratorsGTest, testDynamicBarabasiAlbertGenerator) {
    DynamicGraphSource *gen = new DynamicBarabasiAlbertGenerator(2);

    GraphEventProxy *Gproxy = gen->newGraph();
    Graph *G = Gproxy->G;

    gen->initializeGraph();

    EXPECT_EQ(2u, G->numberOfNodes()) << "initially the generator creates two connected nodes";
    EXPECT_EQ(1u, G->numberOfEdges()) << "initially the generator creates two connected nodes";

    count n = 100;

    gen->generateWhile([&]() { return (G->numberOfNodes() < n); });

    EXPECT_EQ(n, G->numberOfNodes());
    DEBUG("m = ", G->numberOfEdges());

    // resume generator

    gen->generateWhile([&]() { return (G->numberOfNodes() < 2 * n); });
    EXPECT_EQ(2 * n, G->numberOfNodes());

    delete gen;
    delete Gproxy;
    delete G;
}

TEST_F(GeneratorsGTest, viewDynamicBarabasiAlbertGenerator) {
    DynamicGraphSource *gen = new DynamicBarabasiAlbertGenerator(2);
    GraphEventProxy *Gproxy = gen->newGraph();
    Graph *G = Gproxy->G;
    gen->initializeGraph();
    count n = 42;
    gen->generateWhile([&]() { return (G->numberOfNodes() < n); });
    METISGraphWriter writer;
    writer.write(*G, "output/BATest.graph");

    delete gen;
    delete Gproxy;
    delete G;
}

TEST_F(GeneratorsGTest, testStaticPubWebGenerator) {
    Aux::Random::setSeed(42, false);

    count n = 450;
    count numCluster = 9;
    count maxNumNeighbors = 36;
    float rad = 0.075;

    PubWebGenerator gen(n, numCluster, rad, maxNumNeighbors);
    Graph G = gen.generate();
    auto coordinates = gen.moveCoordinates();
    EXPECT_EQ(n, G.numberOfNodes()) << "number of generated nodes";

    // check degree
    G.forNodes([&](node v) { EXPECT_LE(G.degree(v), maxNumNeighbors) << "maximum degree"; });

    // 1-clustering
    ClusteringGenerator clusterGen;
    Partition oneClustering = clusterGen.makeOneClustering(G);
    EXPECT_EQ(G.numberOfNodes(), oneClustering.numberOfElements());

    // output to EPS file
    PostscriptWriter psWriter(true);
    psWriter.write(G, coordinates, oneClustering, "output/pubweb.eps");

    // clustering
    PLM clusterAlgo(G, false, 1.0, "none randomized");
    clusterAlgo.run();
    Partition clustering = clusterAlgo.getPartition();
    EXPECT_EQ(G.numberOfNodes(), clustering.numberOfElements());
    psWriter.write(G, coordinates, clustering, "output/pubweb-clustered-PLM.eps");

    Modularity mod;
    double modVal = mod.getQuality(clustering, G);
    EXPECT_GE(modVal, 0.2) << "modularity of clustering";
    DEBUG("Modularity of clustering: ", modVal);
    DEBUG("Total edge weight: ", G.totalEdgeWeight());
    EXPECT_TRUE(G.checkConsistency());
}

TEST_F(GeneratorsGTest, testDynamicPubWebGenerator) {
    count nSteps = 5;
    count n = 200;
    count numCluster = 10;
    count maxNumNeighbors = 40;
    float rad = 0.08;

    DynamicPubWebGenerator dynGen(n, numCluster, rad, maxNumNeighbors, false);
    Graph G = dynGen.getGraph();
    GraphUpdater gu(G);
    std::vector<GraphEvent> stream;

    for (index i = 1; i <= nSteps; ++i) {
        stream = dynGen.generate(1);
        DEBUG("updating graph");
        gu.update(stream);

        DEBUG("updated graph, new (n, m) = (", G.numberOfNodes(), ", ", G.numberOfEdges(), ")");
        edgeweight tew = G.totalEdgeWeight();
        DEBUG("1/2 graph volume: ", tew);
        EXPECT_GT(tew, 0);

// update coordinates
#ifndef NETWORKIT_RELEASE_LOGGING
        {
            PostscriptWriter psWriter(true);
            std::stringstream ss;
            ss << "output/pubweb-" << std::setw(4) << std::setfill('0') << i << ".eps";
            psWriter.write(G, dynGen.getCoordinates(), ss.str());
        }
#endif
    }
}

/**
 * Testing the dynamic hyperbolic generator with fixed parameters and changing node positions
 */
TEST_F(GeneratorsGTest, testDynamicHyperbolicGeneratorOnMovedNodes) {
    Aux::Random::setSeed(0, false);

    // set up dynamic parameters
    int nSteps = 20;
    const count n = 500;
    const double k = 6;
    const double alpha = 1;
    // const double exp = 2*alpha+1;
    const double T = 0;
    const double R = HyperbolicSpace::getTargetRadius(n, n * k / 2, alpha, T);

    double movedShare = 1;
    double moveDistance = 0.1;

    // set up initial node positions
    vector<double> angles(n, -1);
    vector<double> radii(n, -1);
    HyperbolicSpace::fillPoints(angles, radii, R, alpha);
    DynamicHyperbolicGenerator dynGen(angles, radii, R, alpha, T, movedShare, moveDistance);

    // generate starting graph
    Graph G = HyperbolicGenerator().generate(angles, radii, R);
    count initialEdgeCount = G.numberOfEdges();
    count expected = n * HyperbolicSpace::getExpectedDegree(n, alpha, R) * 0.5;
    EXPECT_NEAR(initialEdgeCount, expected, expected / 5);
    GraphUpdater gu(G);
    std::vector<GraphEvent> stream;

    for (int i = 0; i < nSteps; i++) {
        // move nodes and generate stream of affected edges
        stream = dynGen.generate(1);
        DEBUG("Edges: ", G.numberOfEdges());
        for (auto event : stream) {
            EXPECT_TRUE(event.type == GraphEvent::EDGE_REMOVAL
                        || event.type == GraphEvent::EDGE_ADDITION
                        || event.type == GraphEvent::TIME_STEP);
            if (event.type == GraphEvent::EDGE_REMOVAL) {
                EXPECT_TRUE(G.hasEdge(event.u, event.v));
            }
            // only present nodes can be affected, no new nodes are introduced
            if (event.type != GraphEvent::TIME_STEP) {
                EXPECT_LT(event.u, G.upperNodeIdBound());
            }
        }
        gu.update(stream);
        EXPECT_TRUE(G.checkConsistency());
    }

    // update moved nodes
    angles = getAngles(dynGen);
    radii = getRadii(dynGen);
    Graph comparison = HyperbolicGenerator().generate(angles, radii, R);
    EXPECT_EQ(G.numberOfEdges(), comparison.numberOfEdges());

    // heuristic criterion: Number of edges may change, but should not change much
    EXPECT_NEAR(G.numberOfEdges(), initialEdgeCount, initialEdgeCount / 5);
}

/**
 * creates a series of pictures visualizing the effect of the dynamic hyperbolic generator
 */
TEST_F(GeneratorsGTest, testDynamicHyperbolicVisualization) {
    count n = 300;
    count nSteps = 20;

    const double k = 6;
    const double alpha = 1;
    // const double exp = 2*alpha+1;
    const double T = 0;
    const double R = HyperbolicSpace::getTargetRadius(n, n * k / 2, alpha, T);

    double movedShare = 0.2;
    double moveDistance = 1;
    vector<double> angles(n);
    vector<double> radii(n);

    HyperbolicSpace::fillPoints(angles, radii, R, alpha);

    DynamicHyperbolicGenerator dynGen(angles, radii, R, alpha, T, movedShare, moveDistance);
    Graph G = dynGen.getGraph();

    GraphUpdater gu(G);
    std::vector<GraphEvent> stream;

    for (index i = 0; i < nSteps; i++) {
        stream = dynGen.generate(1);
        DEBUG("Edges: ", G.numberOfEdges());
        for (auto event : stream) {
            EXPECT_TRUE(event.type == GraphEvent::EDGE_REMOVAL
                        || event.type == GraphEvent::EDGE_ADDITION
                        || event.type == GraphEvent::TIME_STEP);
        }
        gu.update(stream);

        auto coordinates = dynGen.getCoordinates();

#ifndef NETWORKIT_RELEASE_LOGGING
        {
            std::stringstream ss;
            ss << "output/hyperbolic-" << std::setw(4) << std::setfill('0') << i << ".eps";

            PostscriptWriter psWriter(true);
            psWriter.write(G, coordinates, ss.str());
        }
#endif
    }
}

TEST_F(GeneratorsGTest, testBarabasiAlbertGeneratorConstructor) {
    // k > nMax
    EXPECT_THROW(BarabasiAlbertGenerator generator(10, 9, 8, true), std::runtime_error);

    // n0 > nMax
    EXPECT_THROW(BarabasiAlbertGenerator generator(5, 9, 10, true), std::runtime_error);

    // n0 = initGraph.numberOfNodes() > nMax
    Graph initGraph(10);
    EXPECT_THROW(BarabasiAlbertGenerator generator(6, 9, initGraph, true), std::runtime_error);

    // initGraph, k > nMax
    initGraph = Graph(6);
    EXPECT_THROW(BarabasiAlbertGenerator generator(10, 9, initGraph, true), std::runtime_error);

    // initGraph does not have consecutive node ids
    initGraph.removeNode(0);
    EXPECT_THROW(BarabasiAlbertGenerator generator(3, 9, initGraph, true), std::runtime_error);
}

TEST_F(GeneratorsGTest, testBarabasiAlbertGeneratorBatagelj) {
    count k = 3;
    count nMax = 1000;
    count n0 = 3;

    BarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0, true);
    Graph G = BarabasiAlbert.generate();

    EXPECT_EQ(nMax, G.numberOfNodes());
    EXPECT_LE(G.numberOfEdges(), nMax * k);
    EXPECT_TRUE(G.checkConsistency());

    Graph initGraph(4);
    initGraph.addEdge(0, 1);
    initGraph.addEdge(2, 1);
    initGraph.addEdge(2, 3);
    initGraph.addEdge(0, 3);
    BarabasiAlbert = BarabasiAlbertGenerator(k, nMax, initGraph, true);
    G = BarabasiAlbert.generate();

    EXPECT_EQ(nMax, G.numberOfNodes());
    EXPECT_LE(G.numberOfEdges(), nMax * k);
    EXPECT_TRUE(G.checkConsistency());
}

TEST_F(GeneratorsGTest, generatetBarabasiAlbertGeneratorGraph) {
    count k = 3;
    count nMax = 1000;
    count n0 = 3;

    BarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0);

    Graph G = BarabasiAlbert.generate();
    GraphIO io;
    io.writeAdjacencyList(G, "output/BarabasiGraph.txt");
}

TEST_F(GeneratorsGTest, testParallelBarabasiAlbertGeneratorDistribution) {
    count k = 3;
    count nMax = 1000;
    count n0 = 3;

    BarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0, false);

    Graph G = BarabasiAlbert.generate();
    G.forNodes([&G, k, n0](node u) {
        if (u < n0)
            return;
        int totalEdges = 0;
        G.forEdgesOf(u, [&totalEdges, u](node v) {
            if (v < u)
                totalEdges++;
        });
        EXPECT_EQ(totalEdges, k);
    });
    PowerlawDegreeSequence sequence(G);
    auto gamma = sequence.getGamma();
    EXPECT_LT(gamma, -2);
    EXPECT_GT(gamma, -3);
}

TEST_F(GeneratorsGTest, testDynamicPathGenerator) {
    count nSteps = 42;
    DynamicPathGenerator gen;
    auto stream = gen.generate(nSteps);
    EXPECT_EQ(stream.size(), nSteps * 3 + 1);
#ifndef NETWORKIT_RELEASE_LOGGING
    for (auto ev : stream) {
        TRACE(ev.toString());
    }
#endif
}

TEST_F(GeneratorsGTest, testErdosRenyiGenerator) {
    count n = 2000;
    double p = 1.5 * (std::log(n) / (double)n);

    ErdosRenyiGenerator generator(n, p);
    Graph G = generator.generate();
    EXPECT_EQ(n, G.numberOfNodes());
    EXPECT_FALSE(G.isEmpty());
    EXPECT_TRUE(G.checkConsistency());

    count nPairs = (n * (n - 1)) / 2;
    count nEdges = G.numberOfEdges();
    EXPECT_GE(nEdges, 0.75 * p * nPairs);
    EXPECT_LE(nEdges, 1.25 * p * nPairs);

    DEBUG("Number of edges with probability ", p, " (actual/expected): ", nEdges, " / ",
          (nPairs * p));
    EXPECT_TRUE(G.checkConsistency());
}

TEST_F(GeneratorsGTest, testRmatGeneratorException) {
    count scale = 9;
    count edgeFactor = 12;
    double a = 0.51;
    double b = 0.12;
    double c = 0.12;
    double d = 0.2;

    EXPECT_THROW(RmatGenerator rmat(scale, edgeFactor, a, b, c, d), std::runtime_error);
}

TEST_F(GeneratorsGTest, testRmatGenerator) {
    count scale = 9;
    count n = (1 << scale);
    count edgeFactor = 12;
    double a = 0.51;
    double b = 0.12;
    double c = 0.12;
    double d = 0.25;

    RmatGenerator rmat(scale, edgeFactor, a, b, c, d);
    Graph G = rmat.generate();

    EXPECT_EQ(G.numberOfNodes(), n);
    EXPECT_LE(G.numberOfEdges(), n * edgeFactor);

    ClusteringCoefficient cc;
    double ccex = cc.exactGlobal(G);
    EXPECT_LE(ccex, 0.4);

    PLM clusterer(G, true);
    clusterer.run();
    Partition zeta = clusterer.getPartition();
    Modularity mod;
    double modVal = mod.getQuality(zeta, G);
    INFO("Modularity of R-MAT graph clustering: ", modVal);
    EXPECT_GE(modVal, 0.0);
    EXPECT_TRUE(G.checkConsistency());
}

TEST_F(GeneratorsGTest, testRmatGeneratorDistribution) {
    constexpr count scale = 2;
    constexpr count n = (1 << scale);
    count edgeFactor = 1;
    double a = 0.51;
    double b = 0.12;
    double c = 0.12;
    double d = 0.25;
    double edgeExpectedProbability[n][n] = {
        {0.0e+00, 0.0e+00, 0.0e+00, 0.0e+00},
        {2.242834007823646e-01, 0.0e+00, 0.0e+00, 0.0e+00},
        {2.242834007823646e-01, 1.0219127659785153e-01, 0.0e+00, 0.0e+00},
        {1.0219127659785153e-01, 1.7352532261978387e-01, 1.7352532261978387e-01, 0.0e+00},
    };

    RmatGenerator rmat(scale, edgeFactor, a, b, c, d, false, 0, true);
    count edgeCount[n][n]{{0}};
    count totalEdges = 0;
    // Now we generate a bunch of graphs and count the edges.
    for (index k = 0; k < 1000; k++) {
        Graph G = rmat.generate();
        G.forEdges([&edgeCount, &totalEdges](node u, node v) {
            edgeCount[u][v] += 1;
            totalEdges += 1;
        });
    }
    for (index i = 0; i < n; ++i) {
        for (index j = 0; j < n; ++j) {
            EXPECT_NEAR((static_cast<double>(edgeCount[i][j]) / static_cast<double>(totalEdges)),
                        edgeExpectedProbability[i][j], 0.01);
        }
    }
}

TEST_F(GeneratorsGTest, testRmatGeneratorReduceNodes) {
    count scale = 9;
    count n = (1 << scale);
    count edgeFactor = 12;
    double a = 0.51;
    double b = 0.12;
    double c = 0.12;
    double d = 0.25;
    int reducedNodes = 4;

    RmatGenerator rmat(scale, edgeFactor, a, b, c, d, false, reducedNodes);
    Graph G = rmat.generate();

    EXPECT_EQ(G.numberOfNodes(), n - reducedNodes);
    EXPECT_LE(G.numberOfEdges(), n * edgeFactor);

    ClusteringCoefficient cc;
    double ccex = cc.exactGlobal(G);
    EXPECT_LE(ccex, 0.4);

    PLM clusterer(G, true);
    clusterer.run();
    Partition zeta = clusterer.getPartition();
    Modularity mod;
    double modVal = mod.getQuality(zeta, G);
    INFO("Modularity of R-MAT graph clustering: ", modVal);
    EXPECT_GE(modVal, 0.0);
    EXPECT_TRUE(G.checkConsistency());
}

TEST_F(GeneratorsGTest, testChungLuGenerator) {
    count n = 400;
    count maxDegree = n / 8;
    std::vector<count> sequence(n);
    count expVolume = 0;
    count actualVolume = 0;

    // fill sequence with random values (this is not power-law, of course!)
    for (index i = 0; i < n; ++i) {
        sequence[i] = rand() % maxDegree;
        expVolume += sequence[i];
    }

    ChungLuGenerator gen(sequence);
    Graph G = gen.generate();
    EXPECT_TRUE(G.checkConsistency());

    EXPECT_EQ(n, G.numberOfNodes());
    G.forNodes([&](node v) { actualVolume += G.degree(v); });

    INFO("expected volume: ", expVolume, ", actual volume: ", actualVolume);
}

TEST_F(GeneratorsGTest, testChungLuGeneratorDegreeConsistency) {
    count n = 1000;
    std::vector<count> vec;
    count maxDegree = n / 8;
    /* Creates a random sequence of weights */
    for (index i = 0; i < n; i++) {
        int grad = Aux::Random::integer(1, maxDegree);
        vec.push_back(grad);
    }
    ChungLuGenerator generator(vec);
    Graph G = generator.generate();
    /* We check to see if the actual degrees of our nodes vary too much from the expected ones.
     * However, we need to sort the expected degrees first, since the algorithm does this as well
     * and the nodes with the highest degrees are added first. */
    Aux::Parallel::sort(vec.begin(), vec.end(), [](count a, count b) { return a > b; });
    /* Check if node degree is more than 50% off from the expected degree of that node. */
    // TODO Should we be looking for something better than a 50% range here?
    G.parallelForNodes([&](node v) { EXPECT_NEAR(G.degree(v), vec[v], (0.5 * maxDegree)); });
}

TEST_F(GeneratorsGTest, testChungLuGeneratorVolumeConsistency) {
    count n = 1000;
    std::vector<count> vec;
    count maxDegree = n / 8;
    count expectedVolume = 0;
    /* Creates a random sequence of weights */
    for (index i = 0; i < n; i++) {
        int grad = Aux::Random::integer(1, maxDegree);
        vec.push_back(grad);
        expectedVolume += grad;
    }
    ChungLuGenerator generator(vec);
    Graph G = generator.generate();
    /* Check if volume is more than 10% off from the expected volume. */
    // TODO Is a 20% offset here sufficient? */
    EXPECT_NEAR(G.numberOfEdges() * 2, expectedVolume, 0.2 * expectedVolume);
}

TEST_F(GeneratorsGTest, testChungLuGeneratorAlamEtAl) {
    count n = 400;
    count maxDegree = n / 8;
    std::vector<count> sequence(n);
    count expVolume = 0;
    count actualVolume = 0;

    // fill sequence with random values (this is not power-law, of course!)
    for (index i = 0; i < n; ++i) {
        sequence[i] = rand() % maxDegree;
        expVolume += sequence[i];
    }

    ChungLuGeneratorAlamEtAl gen(sequence);
    Graph G = gen.generate();
    EXPECT_TRUE(G.checkConsistency());

    EXPECT_EQ(n, G.numberOfNodes());
    G.forNodes([&](node v) { actualVolume += G.degree(v); });

    INFO("expected volume: ", expVolume, ", actual volume: ", actualVolume);
}

TEST_F(GeneratorsGTest, testChungLuGeneratorAlamEtAlDegreeConsistency) {
    count n = 1000;
    std::vector<count> vec;
    count maxDegree = n / 8;
    /* Creates a random sequence of weights */
    for (index i = 0; i < n; i++) {
        int grad = Aux::Random::integer(1, maxDegree);
        vec.push_back(grad);
    }
    ChungLuGeneratorAlamEtAl generator(vec);
    Graph G = generator.generate();
    /* We check to see if the actual degrees of our nodes vary too much from the expected ones.
     * However, we need to sort the expected degrees first, since the algorithm does this as well
     * and the nodes with the highest degrees are added first. */
    Aux::Parallel::sort(vec.begin(), vec.end(), [](count a, count b) { return a < b; });
    /* Check if node degree is more than 50% off from the expected degree of that node. */
    // TODO Should we be looking for something better than a 50% range here?
    G.parallelForNodes([&](node v) { EXPECT_NEAR(G.degree(v), vec[v], (0.5 * maxDegree)); });
}

TEST_F(GeneratorsGTest, testChungLuGeneratorAlamEtAlVolumeConsistency) {
    count n = 1000;
    std::vector<count> vec;
    count maxDegree = n / 8;
    count expectedVolume = 0;
    /* Creates a random sequence of weights */
    for (index i = 0; i < n; i++) {
        int grad = Aux::Random::integer(1, maxDegree);
        vec.push_back(grad);
        expectedVolume += grad;
    }
    ChungLuGeneratorAlamEtAl generator(vec);
    Graph G = generator.generate();
    /* Check if volume is more than 10% off from the expected volume. */
    // TODO Is a 20% offset here sufficient? */
    EXPECT_NEAR(G.numberOfEdges() * 2, expectedVolume, 0.2 * expectedVolume);
}

TEST_F(GeneratorsGTest, testHavelHakimiGeneratorOnRandomSequence) {
    count n = 400;
    std::vector<count> sequence(n);
    bool realizable = false;

    std::mt19937_64 prng(1);
    std::uniform_int_distribution<count> deg_distr(0, n / 10 - 1);
    do {
        // fill sequence with random values (this is not power-law, of course!)
        std::generate(sequence.begin(), sequence.end(), [&] { return deg_distr(prng); });

        // check if sequence is realizable
        HavelHakimiGenerator hhgen(sequence);
        realizable = hhgen.isRealizable();

        if (realizable) {
            Graph G = hhgen.generate();
            EXPECT_TRUE(G.checkConsistency());
            count volume = std::accumulate(sequence.begin(), sequence.end(), 0);
            EXPECT_EQ(volume, 2 * G.numberOfEdges());
        }
    } while (!realizable);
}

TEST_F(GeneratorsGTest, testHavelHakimiGeneratorOnRealSequence) {
    METISGraphReader reader;
    auto graphs = {"input/jazz.graph", "input/lesmis.graph"};

    for (auto path : graphs) {
        Graph G = reader.read(path);
        count n = G.numberOfNodes();
        std::vector<count> sequence(n);
        G.forNodes([&](node u) { sequence[u] = G.degree(u); });

        HavelHakimiGenerator hhgen(sequence);
        Graph G2 = hhgen.generate();
        EXPECT_TRUE(G.checkConsistency());

        count volume = std::accumulate(sequence.begin(), sequence.end(), 0);
        EXPECT_EQ(volume, 2 * G2.numberOfEdges());

        if (volume < 50000) {
            std::vector<count> testSequence(n);
            G2.forNodes([&](node u) { testSequence[u] = G2.degree(u); });

            for (index i = 0; i < n; ++i) {
                EXPECT_EQ(sequence[i], testSequence[i]);
            }
        }
    }
}

TEST_F(GeneratorsGTest, testHavelHakimiGeneratorOnUnrealizableSequence) {
    std::vector<count> seq = {20, 10, 2, 2, 2, 2, 2, 2, 2, 2, 2};

    HavelHakimiGenerator hhgen(seq);
    EXPECT_THROW(hhgen.generate(), std::runtime_error);

    hhgen = HavelHakimiGenerator(seq, true);
    Graph G = hhgen.generate();

    G.forNodes([&](node u) { EXPECT_EQ(std::min<count>(seq[u], 10), G.degree(u)); });
}

TEST_F(GeneratorsGTest, testDynamicForestFireGenerator) {
    Graph G1(0);
    GraphUpdater gu1(G1);
    std::vector<GraphEvent> stream;
    DynamicForestFireGenerator ffg1(0.0, false);
    stream = ffg1.generate(10);
    gu1.update(stream);
    EXPECT_TRUE(G1.checkConsistency());
    EXPECT_EQ(10u, G1.numberOfNodes());
    G1.forNodes([&](node u) {
        count c = 0;
        G1.forNeighborsOf(u, [&](node v) {
            if (v < u) {
                c += 1;
            }
        });
        if (u == 0) {
            EXPECT_EQ(0u, c);
        } else {
            EXPECT_EQ(1u, c);
        }
    });

    Graph G2(0);
    GraphUpdater gu2(G2);
    DynamicForestFireGenerator ffg2(1.0, true, 1.0);
    stream = ffg2.generate(10);
    gu2.update(stream);
    EXPECT_TRUE(G2.checkConsistency());
    EXPECT_EQ(10u, G2.numberOfNodes());
    G2.forNodePairs([&](node u, node v) {
        if (v < u) {
            EXPECT_TRUE(G2.hasEdge(u, v));
        }
    });
    stream = ffg2.generate(10);
    gu2.update(stream);
    EXPECT_EQ(20u, G2.numberOfNodes());
}

TEST_F(GeneratorsGTest, testRegularRingLatticeGenerator) {
    int n0 = 10;
    int neighbors = 2;
    auto testRingLattice = [&](Graph G) {
        EXPECT_EQ(n0, (int)G.numberOfNodes());
        EXPECT_EQ(n0 * neighbors, (int)G.numberOfEdges());
        G.forNodePairs([&](node u, node v) {
            int diff = std::abs((int)u - (int)v);
            if (u != v && (diff <= neighbors || diff >= n0 - neighbors)) {
                EXPECT_TRUE(G.hasEdge(u, v));
            } else {
                EXPECT_FALSE(G.hasEdge(u, v));
            }
        });
    };

    RegularRingLatticeGenerator rrlg = RegularRingLatticeGenerator(n0, neighbors);
    testRingLattice(rrlg.generate());
}

TEST_F(GeneratorsGTest, testWattsStrogatzGenerator) {
    int n0 = 10;
    int neighbors = 2;
    auto testRingLattice = [&](Graph G) {
        G.forNodePairs([&](node u, node v) {
            int diff = std::abs((int)u - (int)v);
            if (u != v && (diff <= neighbors || diff >= n0 - neighbors)) {
                EXPECT_TRUE(G.hasEdge(u, v));
            } else {
                EXPECT_FALSE(G.hasEdge(u, v));
            }
        });
    };

    WattsStrogatzGenerator wsg1 = WattsStrogatzGenerator(n0, neighbors, 0.0);
    testRingLattice(wsg1.generate());

    WattsStrogatzGenerator wsg2 = WattsStrogatzGenerator(n0, neighbors, 0.3);
    Graph G = wsg2.generate();
    EXPECT_TRUE(G.checkConsistency());
    EXPECT_EQ(n0, (int)G.numberOfNodes());
    EXPECT_EQ(n0 * neighbors, (int)G.numberOfEdges());
}

TEST_F(GeneratorsGTest, testWattsStrogatzGeneratorBigKs) {
    constexpr count nodes = 10;
    constexpr count neighbors = 4;
    constexpr double p = 0.5;
    for (int seed : {1, 2, 3, 4, 5}) {
        Aux::Random::setSeed(seed, false);
        const auto G = WattsStrogatzGenerator(nodes, neighbors, p).generate();
        EXPECT_TRUE(G.checkConsistency());
        EXPECT_EQ(nodes, G.numberOfNodes());
        EXPECT_EQ(nodes * neighbors, G.numberOfEdges());
    }
}

TEST_F(GeneratorsGTest, testDorogovtsevMendesGenerator) {
    int n0 = 20;
    DorogovtsevMendesGenerator dmg = DorogovtsevMendesGenerator(n0);
    Graph G = dmg.generate();

    EXPECT_EQ(n0, (int)G.numberOfNodes());
    EXPECT_EQ(2 * n0 - 3, (int)G.numberOfEdges());
    G.forNodes([&](node u) {
        count c = 0;
        G.forNeighborsOf(u, [&](node v) {
            if (v < u) {
                c += 1;
            }
        });
        if (u <= 2) {
            EXPECT_EQ(u, c);
        } else {
            EXPECT_EQ(2u, c);
        }
    });
    EXPECT_TRUE(G.checkConsistency());
}

TEST_F(GeneratorsGTest, testDynamicDorogovtsevMendesGenerator) {
    count n0 = 20;
    DynamicDorogovtsevMendesGenerator ddmg = DynamicDorogovtsevMendesGenerator();
    Graph G(0);
    GraphUpdater gu(G);
    std::vector<GraphEvent> stream;
    stream = ddmg.generate(n0 - 3);
    gu.update(stream);

    EXPECT_EQ(n0, G.numberOfNodes());
    EXPECT_EQ(2 * n0 - 3, G.numberOfEdges());
    G.forNodes([&](node u) {
        count c = 0;
        G.forNeighborsOf(u, [&](node v) {
            if (v < u) {
                c += 1;
            }
        });
        if (u <= 2) {
            EXPECT_EQ(u, c);
        } else {
            EXPECT_EQ(2u, c);
        }
    });
}

TEST_F(GeneratorsGTest, testStaticDegreeSequenceGenerator) {
    auto test_known = [](const std::vector<count> &seq, bool result) {
        HavelHakimiGenerator gen(seq, true);
        ASSERT_EQ(gen.isRealizable(), result);
    };

    // some contrived sequence
    test_known({1}, false);
    test_known({1, 1}, true);
    test_known({2, 2, 2}, true);
    test_known({1, 2, 1}, true);
    test_known({1, 3, 1}, false);
    test_known({1, 1, 3, 1}, true);

    // some random sequences
    std::mt19937_64 prng(1);
    std::uniform_int_distribution<node> distr_num_nodes(1, 50);

    unsigned numRealized = 0;
    for (int iter = 0; iter < 100; ++iter) {
        const auto n = distr_num_nodes(prng);
        std::vector<count> seq(n);
        std::generate(seq.begin(), seq.end(), [&] {
            return std::uniform_int_distribution<count>{0, n - 1}(prng);
        });

        HavelHakimiGenerator gen(seq, true);
        const auto isRealizable = gen.isRealizable();

        const auto G = gen.generate();
        const auto nodes = G.nodeRange();
        ASSERT_EQ(G.numberOfNodes(), n);

        const auto didRealize =
            std::all_of(nodes.begin(), nodes.end(), [&](node u) { return G.degree(u) == seq[u]; });

        ASSERT_EQ(isRealizable, didRealize);
        numRealized += didRealize;
    }

    ASSERT_GT(numRealized, 10);
    ASSERT_LT(numRealized, 90);
}

TEST_F(GeneratorsGTest, testStochasticBlockmodel) {
    count n = 10;
    count nBlocks = 2;
    std::vector<index> membership = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
    std::vector<std::vector<double>> affinity = {{1.0, 0.0}, {0.0, 1.0}};
    StochasticBlockmodel sbm(n, nBlocks, membership, affinity);
    Graph G = sbm.generate();

    EXPECT_EQ(n, G.numberOfNodes());
    EXPECT_EQ(20u, G.numberOfEdges());
}

/**
 * Test whether points generated in hyperbolic space fulfill basic constraints
 */
TEST_F(GeneratorsGTest, testHyperbolicPointGeneration) {
    count n = 1000;
    double stretch = Aux::Random::real(0.5, 1.5);
    double alpha = Aux::Random::real(0.5, 1.5);
    double R = HyperbolicSpace::hyperbolicAreaToRadius(n) * stretch;
    vector<double> angles(n, -1);
    vector<double> radii(n, -1);
    HyperbolicSpace::fillPoints(angles, radii, R, alpha);
    for (index i = 0; i < n; i++) {
        EXPECT_GE(angles[i], 0);
        EXPECT_LT(angles[i], 2 * PI);
        EXPECT_GE(radii[i], 0);
        EXPECT_LE(radii[i], R);
    }
}

/**
 * Test whether the number edges generated by the hyperbolic generator agree at least roughly with
 * theory
 */
TEST_F(GeneratorsGTest, testHyperbolicGenerator) {
    Aux::Random::setSeed(0, false);
    count n = 5000;
    double k = 16;
    count m = k * n / 2;
    HyperbolicGenerator gen(n, k, 7);
    Graph G = gen.generate();
    EXPECT_EQ(G.numberOfNodes(), n);
    EXPECT_TRUE(G.checkConsistency());
    EXPECT_NEAR(G.numberOfEdges(), m, m / 5);
}

/**
 * Check consistency of graphs generated by the hyperbolic generator
 */
TEST_F(GeneratorsGTest, testHyperbolicGeneratorConsistency) {
    Aux::Random::setSeed(0, false);
    count n = 5000;
    double k = 6;
    count m = n * k / 2;
    HyperbolicGenerator gen(n, k);
    Graph G = gen.generate();
    EXPECT_NEAR(G.numberOfEdges(), m, m / 5);
    ASSERT_TRUE(G.checkConsistency());
}

TEST_F(GeneratorsGTest, testHyperbolicGeneratorMechanicGraphs) {
    Aux::Random::setSeed(0, false);
    count n = 2000;
    double k = 6;
    count m = n * k / 2;
    HyperbolicGenerator gen(n, k, 3, 0.14);
    Graph G = gen.generate();
    EXPECT_NEAR(G.numberOfEdges(), m, m / 10);
    ASSERT_TRUE(G.checkConsistency());
}

TEST_F(GeneratorsGTest, testConfigurationModelGeneratorOnRealSequence) {
    METISGraphReader reader;
    auto graphs = {"input/jazz.graph", "input/lesmis.graph"};

    for (auto path : graphs) {
        Graph G = reader.read(path);
        count n = G.numberOfNodes();
        std::vector<count> sequence(n);
        G.forNodes([&](node u) { sequence[u] = G.degree(u); });

        bool skipTest = false;
        EdgeSwitchingMarkovChainGenerator gen(sequence, skipTest);
        Graph G2 = gen.generate();

        count volume = std::accumulate(sequence.begin(), sequence.end(), 0);
        EXPECT_EQ(volume, 2 * G2.numberOfEdges());

        if (volume < 50000) {
            std::vector<count> testSequence(n);
            G2.forNodes([&](node u) { testSequence[u] = G2.degree(u); });
            Aux::Parallel::sort(testSequence.begin(), testSequence.end(), std::greater<count>());
            Aux::Parallel::sort(sequence.begin(), sequence.end(), std::greater<count>());

            for (index i = 0; i < n; ++i) {
                EXPECT_EQ(sequence[i], testSequence[i]);
            }
        }
    }
}

TEST_F(GeneratorsGTest, debugHyperbolicHighTemperatureGraphs) {
    count n = 10000;
    double k = 10;
    double gamma = 3;
    count m = n * k / 2;
    for (double T = 0; T < 10; T += 0.1) {
        if (std::abs(T - 1) < 0.00001)
            continue;
        HyperbolicGenerator gen(n, k, gamma, T);
        Graph G = gen.generate();
        EXPECT_NEAR(G.numberOfEdges(), m, m / 10);
    }
}

TEST_F(GeneratorsGTest, debugGiganticCollectionOfHyperbolicTemperatureGraphs) {
    for (index i = 0; i < 30; i++) {
        count n = 10000;
        double k = 10;
        double T = 0.1;
        count m = n * k / 2;
        HyperbolicGenerator gen(n, k, 3, T);
        Graph G = gen.generate();
        EXPECT_NEAR(G.numberOfEdges(), m, m / 10);
    }
}

TEST_F(GeneratorsGTest, debugGiganticCollectionOfHyperbolicUnitDiskGraphs) {
    count n = 1000000;
    double k = 1;
    for (index i = 0; i < 7; i++) {
        count m = n * k / 2;
        HyperbolicGenerator gen(n, k, 7);
        Graph G = gen.generate();
        EXPECT_NEAR(G.numberOfEdges(), m, m / 5);
        EXPECT_TRUE(G.checkConsistency());
        k *= 2;
    }
}

TEST_F(GeneratorsGTest, testLFRGenerator) {
    Aux::Random::setSeed(42, false);
    count n = 500;
    LFRGenerator gen(n);
    gen.generatePowerlawDegreeSequence(20, 50, -2);
    gen.generatePowerlawCommunitySizeSequence(10, 50, -1);
    gen.setMu(0.5);
    gen.run();
    Graph G1 = gen.getMoveGraph();
    gen.run(); // should rewire the edges but nothing else
    Graph G2 = gen.getMoveGraph();
    EXPECT_EQ(n, G1.numberOfNodes());
    EXPECT_EQ(n, G2.numberOfNodes());
    EXPECT_EQ(G1.numberOfEdges(), G2.numberOfEdges());
}

TEST_F(GeneratorsGTest, testLFRGeneratorImpossibleSequence) {
    LFRGenerator gen(100);
    gen.generatePowerlawDegreeSequence(10, 11, -2);
    EXPECT_ANY_THROW(gen.generatePowerlawCommunitySizeSequence(9, 8, -3));
    gen.setMu(0.5);
    EXPECT_THROW(gen.run(), std::runtime_error);
    EXPECT_THROW(gen.getMoveGraph(), std::runtime_error);
}

TEST_F(GeneratorsGTest, testLFRGeneratorWithRealData) {
    const auto degreeSequence = readVector<count>("input/testLFR/testLFRDegSeq.txt"),
               partition = readVector<count>("input/testLFR/testLFRPartition.txt");

    const auto mu = readVector<double>("input/testLFR/testLFRMu.txt");
    Partition C(partition.size());
    C.setUpperBound(20);
    for (node u = 0; u < partition.size(); ++u) {
        C[u] = partition[u];
    }
    LFRGenerator gen(degreeSequence.size());
    gen.setDegreeSequence(degreeSequence);
    gen.setPartition(C);
    gen.setMu(mu);
    gen.run();
    Graph G = gen.getGraph();
    G.parallelForNodes([&](node u) { EXPECT_EQ(G.degree(u), degreeSequence[u]); });
    EXPECT_EQ(C.numberOfSubsets(), gen.getPartition().numberOfSubsets());
}

TEST_F(GeneratorsGTest, testLFRGeneratorWithBinomialDistribution) {
    Aux::Random::setSeed(42, false);
    count n = 500;
    LFRGenerator gen(n);
    gen.generatePowerlawDegreeSequence(20, 50, -2);
    gen.generatePowerlawCommunitySizeSequence(10, 50, -1);
    gen.setMuWithBinomialDistribution(0.5);
    gen.run();
    Graph G1 = gen.getMoveGraph();
    gen.run(); // should rewire the edges but nothing else
    Graph G2 = gen.getMoveGraph();
    EXPECT_EQ(n, G1.numberOfNodes());
    EXPECT_EQ(n, G2.numberOfNodes());
    EXPECT_EQ(G1.numberOfEdges(), G2.numberOfEdges());
}

TEST_F(GeneratorsGTest, testMocnikGenerator) {
    count dim = 3;
    count n = 10000;
    double k = 2.6;

    MocnikGenerator Mocnik(dim, n, k);
    Graph G(0);
    EXPECT_TRUE(G.isEmpty());
    G = Mocnik.generate();
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(G.numberOfNodes(), n);
    EXPECT_NEAR(G.numberOfEdges() * 1. / G.numberOfNodes(), std::pow(k, dim), 20000);
}

TEST_F(GeneratorsGTest, testMocnikGeneratorBasic) {
    count dim = 3;
    count n = 5000;
    double k = 2.6;

    MocnikGeneratorBasic Mocnik(dim, n, k);
    Graph G(0);
    EXPECT_TRUE(G.isEmpty());
    G = Mocnik.generate();
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(G.numberOfNodes(), n);
    EXPECT_NEAR(G.numberOfEdges() * 1. / G.numberOfNodes(), std::pow(k, dim), 10000);
}

TEST_F(GeneratorsGTest, testPowerLawDegreeSequenceFromDegreeSequence) {
    std::vector<double> inputVector = {2, 4, 3, 2, 1, 6};
    double sum = 0.0;
    for (auto it = inputVector.begin(); it != inputVector.end(); ++it) {
        sum += *it;
    }
    auto avg = sum / inputVector.size();

    PowerlawDegreeSequence PLDS(inputVector);
    PLDS.setGammaFromAverageDegree(avg);
    PLDS.run();
    EXPECT_NEAR(PLDS.getExpectedAverageDegree(), avg, 1.0);
    EXPECT_EQ(PLDS.getMinimumDegree(), 1);
    EXPECT_EQ(PLDS.getMaximumDegree(), 6);
}

TEST_F(GeneratorsGTest, testPowerLawDegreeSequenceFromGraph) {
    METISGraphReader reader;
    auto graphpath = "input/jazz.graph";
    Graph G = reader.read(graphpath);

    PowerlawDegreeSequence PLDS(G);
    PLDS.run();
    EXPECT_NEAR(PLDS.getExpectedAverageDegree(), 19, 1);
    EXPECT_EQ(PLDS.getMinimumDegree(), 1);
    EXPECT_EQ(PLDS.getMaximumDegree(), 100);
}

} /* namespace NetworKit */

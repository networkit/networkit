/*
 * FiniteEmbeddingTest.cpp
 *
 *  Created on: 16.10.2020
 *      Author: Klaus Ahrens
 */

#include <iomanip>
#include <iostream>

#include <gtest/gtest.h>

/*
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/centrality/ApproxBetweenness.hpp>
#include <networkit/centrality/ApproxCloseness.hpp>
#include <networkit/centrality/ApproxSpanningEdge.hpp>
#include <networkit/centrality/ApproxGroupBetweenness.hpp>
#include <networkit/centrality/Betweenness.hpp>
#include <networkit/centrality/Closeness.hpp>
#include <networkit/centrality/CoreDecomposition.hpp>
#include <networkit/centrality/DynApproxBetweenness.hpp>
#include <networkit/centrality/DynKatzCentrality.hpp>
#include <networkit/centrality/DynTopHarmonicCloseness.hpp>
#include <networkit/centrality/EigenvectorCentrality.hpp>
#include <networkit/centrality/EstimateBetweenness.hpp>
#include <networkit/centrality/GedWalk.hpp>
#include <networkit/centrality/GroupCloseness.hpp>
#include <networkit/centrality/GroupDegree.hpp>
#include <networkit/centrality/HarmonicCloseness.hpp>
#include <networkit/centrality/KPathCentrality.hpp>
#include <networkit/centrality/KadabraBetweenness.hpp>
#include <networkit/centrality/KatzCentrality.hpp>
#include <networkit/centrality/LaplacianCentrality.hpp>
#include <networkit/centrality/LocalClusteringCoefficient.hpp>
#include <networkit/centrality/PageRank.hpp>
#include <networkit/centrality/PermanenceCentrality.hpp>
#include <networkit/centrality/SpanningEdgeCentrality.hpp>
#include <networkit/centrality/TopCloseness.hpp>
#include <networkit/centrality/TopHarmonicCloseness.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/generators/DorogovtsevMendesGenerator.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/SNAPGraphReader.hpp>
#include <networkit/structures/Cover.hpp>
#include <networkit/structures/Partition.hpp>
*/

#include <networkit/io/METISGraphReader.hpp>
#include <networkit/embedding/Node2Vec.hpp>

namespace NetworKit {

class FiniteEmbeddingTest : public testing::Test{};

using Embeddings = std::vector<std::vector<float>>;

bool allFinite(const Embeddings& e) {
    for (auto& emb: e) {
        for(auto f: emb) {
            if (std::isinf(f)) {
                return false;
            }
        }
    }
    return true;
}

TEST_F(FiniteEmbeddingTest, testFiniteEmbeddingOnSmallGraph) {

    NetworKit::METISGraphReader reader;

    NetworKit::Graph graph = reader.read("input/karate.graph");
    
    auto algo = NetworKit::Node2Vec(graph, .3, 1.3);

    algo.run();
    
    auto features = algo.getFeatures();

    EXPECT_TRUE(allFinite(features));
}

} // namespace NetworKit

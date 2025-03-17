/*
 *  KHop.cpp
 *
 *
 *  Created on: 23.06.2024
 *      Author: Alexander K. Ziebs
 *
 */

#ifndef NETWORKIT_EMBEDDING_K_HOP_HPP_
#define NETWORKIT_EMBEDDING_K_HOP_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup embedding
 */
class KHop final : public Algorithm {

public:
    enum class KHopMode : unsigned char {
        STRICT = 0,
        DEFAULT = 1,
    };

    /**
     * Generates k-hop Graph and probabilities as weights for walks.
     *
     * @param G   The Graph.
     * @param K   Number of hops.
     * @param L   Walk length.
     * @param N   Walk count.
     * @param D   Dimension of embedding vectors.
     * @param M   STRICT (0) for k-hops exactly k, DEFAULT (1) for additional edges less than k.
     * @param winSize    Window Size for walks.
     * @param iterations Iterations for walks.
     */
    KHop(const Graph &G, size_t K = 2, count L = 80, count N = 10, count D = 128,
         KHopMode M = KHopMode::DEFAULT, count winSize = 8, count iterations = 20);

    ~KHop() override = default;

    /**
     This method computes all node embeddings.
    */
    void run() override;

    /**
     This method returns a vector that contains feature vectors for all nodes
    */
    const std::vector<std::vector<float>> &getFeatures() const;

    /**
     * Returns the computed k-hop graph.
     *
     * @return Graph
     */
    Graph getKHopGraph() const { return kHopGraph; }

private:
    // The graph
    const Graph &G;
    // Neighbourhood size K
    size_t K;
    // Walk length L
    count L;
    // Walk count N
    count N;
    // Feature Dimension D
    count D;
    // Mode for k-hop graph generation
    KHopMode M;
    // Windowsize walk
    count winSize;
    // Iterations walk
    count iterations;

    std::vector<std::vector<float>> features;
    Graph kHopGraph;

    void computeKHopGraph();
    void generateEdgeWeights();
};

} /* namespace NetworKit */

#endif // NETWORKIT_EMBEDDING_K_HOP_HPP_

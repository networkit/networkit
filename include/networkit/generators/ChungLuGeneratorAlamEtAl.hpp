/*
 * ChungLu.hpp
 *
 *  Created on: Dec 23, 2013
 *      Author: Henning
 */

#ifndef NETWORKIT_GENERATORS_CHUNG_LU_GENERATOR_ALAM_ET_AL_HPP_
#define NETWORKIT_GENERATORS_CHUNG_LU_GENERATOR_ALAM_ET_AL_HPP_

#include <networkit/generators/StaticDegreeSequenceGenerator.hpp>
#include <networkit/generators/StaticGraphGenerator.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 * Given an arbitrary degree sequence, the Chung-Lu generative model
 * will produce a random graph with the same expected degree sequence.
 *
 * see Chung, Lu: The average distances in random graphs with given expected degrees
 * and Chung, Lu: Connected Components in Random Graphs with Given Expected Degree Sequences.
 * Aiello, Chung, Lu: A Random Graph Model for Massive Graphs describes a different generative model
 * which is basically asymptotically equivalent but produces multi-graphs.
 *
 * This follows the implementation of Joel Miller and Aric Hagberg's
 * "Efficient Generation of Networks with Given Expected Degrees" (2011)
 * http://aric.hagberg.org/papers/miller-2011-efficient.pdf .
 * It gives a complexity of O(n+m) as opposed to quadratic.
 *
 */

class ChungLuGeneratorAlamEtAl final : StaticGraphGenerator {
    struct VertexGroup {
        count degrees;
        count vertexCount;
        count startIndex;
    };
    std::vector<VertexGroup> groups;
    count sum_deg;
    count n;
    bool parallel;

    template <typename F>
    void edgeSkipping(std::mt19937_64 &generator, F &&addEdge, index i, index j, double p,
                      index end);
    Graph generateSequential();
    Graph generateParallel();

public:
    ChungLuGeneratorAlamEtAl(const std::vector<count> &degreeSequence, bool parallel = false);

    /**
     * Generates graph with expected degree sequence seq.
     */
    Graph generate() override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_GENERATORS_CHUNG_LU_GENERATOR_ALAM_ET_AL_HPP_
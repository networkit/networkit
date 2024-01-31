/*
 * ChungLuGeneratorAlamEtAl.hpp
 *
 *  Created on: Dec 23, 2013
 *      Author: Henning
 *      Contributors: Hoske/Weisbarth/Hering
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
 * For reference of the concept of the original model see Chung, Lu: The average distances in random
 * graphs with given expected degrees and Chung, Lu: Connected Components in Random Graphs with
 * Given Expected Degree Sequences. Aiello, Chung, Lu: A Random Graph Model for Massive Graphs
 * describes a different generative model which is basically asymptotically equivalent but produces
 * multi-graphs.
 *
 * This implementation follows a newer parallelizable algorithm by M. Alam, M. Khan, A. Vullikanti
 * and M. Marathe, "An Efficient and Scalable Algorithmic Method for Generating Large-Scale Random
 * Graphs". https://ieeexplore.ieee.org/document/7877110 . The time complexity is O(m/p + delta +
 * p), where m is the number of edges, p are the number of threads, and delta is the number of
 * distinct node degrees from the input sequence. Furthermore the space complexity is O(delta) for
 * each of the p threads.
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

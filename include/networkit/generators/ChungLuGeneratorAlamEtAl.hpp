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
 * The Chung-Lu generative model produces a random graph for a given expected degree sequence.
 * For reference, see the original model proposed by Chung and Lu in [1] and [2].
 * This implementation follows a newer parallelizable algorithm [3] by M. Alam, M. Khan et al.:
 * The time complexity is O(m/p + delta + p), where m is the number of edges, p are the number of
 * threads, and delta is the number of distinct node degrees from the input sequence. Furthermore
 * the space complexity is O(delta) for each of the p threads.
 *
 * [1] F. Chung, L. Lu: [The Average Distance in a Random Graph with Given Expected
 * Degree](https://www.researchgate.net/publication/11004006_The_Average_Distance_in_a_Random_Graph_with_Given_Expected_Degree)
 *
 * [2] F. Chung, L. Lu: [Connected Components in Random Graphs with Given Expected Degree
 * Sequences](https://link.springer.com/article/10.1007/PL00012580)
 *
 * [3] M. Alam, M. Khan, et al.: [An Efficient and Scalable Algorithmic Method for Generating
 * Large-Scale Random Graphs](https://ieeexplore.ieee.org/document/7877110)
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

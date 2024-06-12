/*
 * RmatGenerator.hpp
 *
 *  Created on: 18.03.2014
 *      Author: Henning
 */

#ifndef NETWORKIT_GENERATORS_RMAT_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_RMAT_GENERATOR_HPP_

#include <stdint.h>
#include <networkit/generators/StaticGraphGenerator.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 * Generates static R-MAT graphs. R-MAT (recursive matrix) graphs are
 * random graphs with n=2^scale nodes and m=n*edgeFactor edges.
 * More details at http://www.graph500.org or in the original paper:
 * Deepayan Chakrabarti, Yiping Zhan, Christos Faloutsos:
 * R-MAT: A Recursive Model for Graph Mining. SDM 2004: 442-446.
 */
class RmatGenerator final : public StaticGraphGenerator<Graph> {
    count scale; ///< n = 2^scale
    count edgeFactor;
    double defaultEdgeWeight;
    bool weighted;
    count reduceNodes;
    bool discardSelfLoops;

    // Data for the alias table:
    std::vector<std::pair<uint32_t, uint32_t>> bits;
    std::vector<uint8_t> numberOfBits;
    std::vector<uint32_t> coinFlipProbability;
    std::vector<uint32_t> coinFlipReplacement;
    uint32_t mask;
    std::pair<uint32_t, uint32_t> curBits{0, 0};
    uint32_t remainingBits = 0;

public:
    /**
     * @param[in] scale Number of nodes = 2^scale
     * @param[in] edgeFactor Number of edges = number of nodes * edgeFactor
     * @param[in] a Probability for quadrant upper left
     * @param[in] b Probability for quadrant upper right
     * @param[in] c Probability for quadrant lower left
     * @param[in] d Probability for quadrant lower right
     * @param[in] weighted  result graph weighted?
     * @param[in] reduceNodes  number of random nodes to delete to achieve a given node count
     * @param[in] discardSelfLoops ignore self loops
     */
    RmatGenerator(count scale, count edgeFactor, double a, double b, double c, double d,
                  bool weighted = false, count reduceNodes = 0, bool discardSelfLoops = true);

    /**
     * @return Graph to be generated according to parameters specified in constructor.
     */
    Graph generate() override;

private:
    std::pair<node, node> sampleEdge(uint8_t bits);
};

} /* namespace NetworKit */
#endif // NETWORKIT_GENERATORS_RMAT_GENERATOR_HPP_

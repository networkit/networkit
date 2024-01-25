/*
 *  ComplexPaths.hpp
 *
 *
 *  Created on:16.06.2023
 *      Author: Klaus Ahrens
 *              <ahrens@informatik.hu-berlin.de>
 *
 *  adapted and reimplemented from
 *
 *            https://github.com/drguilbe/complexpaths.git
 *
 *  see [ Guilbeault, D., Centola, D. Topological measures for
 *        identifying and predicting the spread of complex contagions.
 *        Nat Commun 12, 4430 (2021).
 *        https://doi.org/10.1038/s41467-021-24704-6 ]
 *
 */

#ifndef NETWORKIT_CENTRALITY_COMPLEX_PATHS_HPP_
#define NETWORKIT_CENTRALITY_COMPLEX_PATHS_HPP_

#include <random>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 * computes
 *      in Mode::singleNode :
 *	complex path graphs (with inner connection degree above threshold)
 *	from a start node in Mode::singleNode
 *	OR
 *      in Mode::allNodes :
 *	complex path lengths (maximal distance in complex path graphs for
 *	starting nodes u) for all nodes u
 *
 */

class ComplexPathAlgorithm : public Algorithm {
public:
    /**
     * operation modes
     */
    enum Mode { singleNode, allNodes };

    /**
     * Constructs the ComplexPathAlgorithm class for the given Graph @a G.
     * depending on the @a mode the algorithm the algorithm
     *   in @a mode <code>Mode::singleNode</code> starting from @a start
     *      constructs a subgraph of @a G in which all inner nodes have
     *      more than @a threshold neighbors
     *      call <code>getComplexGraph</code> after run to get it
     *   in @a mode <code>Mode::allNodes</code> constructs these subgraphs for all nodes
     *      and calculates their longest paths from their start nodes
     *      these length can be absolute or normalized to [0..1] by calling
     *      normalize before or after <code>run</code>
     *      call <code>getPLci</code> after run to get these lengths
     *
     * @param G The graph.
     * @param threshold number of neighbors needed
     * @param mode as explained above
     * @param start start node for Mode::singleNode
     */
    ComplexPathAlgorithm(const Graph &G, count threshold = 3, Mode mode = Mode::allNodes,
                         node start = none);

    /**
     * normalize path lengths
     */
    void normalize();

    void run() override;

    /**
     * [normalized] results after running in <code>Mode::allNodes</code>
     */
    std::vector<double> getPLci();

    /**
     * resulting graph after running in <code>Mode::singleNode</code>
     */
    Graph getComplexGraph();

    /**
     * nodes in the resulting graph after running in <code>Mode::singleNode</code>
     * which are connected to start by at least threshold paths
     */
    std::vector<node> getAdopters();

private:
    const Graph *inputGraph;
    Graph complexPathGraph;
    std::vector<double> complexPathsLengths;
    std::vector<node> adopters;
    const Mode mode;
    const node start;
    const count threshold;
    bool normPaths;

    Graph complexPathsGraph(node seed, count threshold, std::vector<node> *adopters);
    std::vector<double> complexPathLength(count t);
    std::vector<node> generateSeeds(node seed, const Graph &g, count threshold);
};

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_COMPLEX_PATHS_HPP_

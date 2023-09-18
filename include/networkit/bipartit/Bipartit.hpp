/*
 * Bipartit.hpp
 *
 * Created on: 18.09.2023
 *     Author: Michael Kaibel
 */

#ifndef NETWORKIT_BIPARTIT_BIPARTIT_HPP_
#define NETWORKIT_BIPARTIT_BIPARTIT_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

class Bipartit final : public Algorithm {
public:
    /**
     * Creates Bipartit class for graph @G
     *
     * @param G the graph
     */
    Bipartit(const Graph &G);

    /**
     * Calculates if graph @G is bipartit. If yes generates bipartition, if no finds odd cycle
     */
    void run() override;

    /**
     * Check if the graph is bipartit
     *
     * @return true iff the graph is bipartit
     */
    bool isBipartit();

    /**
     * Get a bipartition iff the graph is bipartit
     *
     * @return a bipartition iff the graph is bipartit, throws error otherwise
     */
    const Partition &getPartition();

    /**
     * Get an odd cycle iff one exists in the graph
     *
     * @return on odd cycle iff the graph has one, throws error otherwise
     */
    const std::vector<node> &getOddCycle();

protected:
    const Graph *G;

    Partition partition;

    std::vector<node> oddCircle;

    bool bipartit = false;

    void findOddCircle(std::vector<node> &parent, node v, node w);
};

} // NetworKit

#endif // NETWORKIT_BIPARTIT_BIPARTIT_HPP_

/*
 * Bipartite.hpp
 *
 * Created on: 18.09.2023
 *     Author: Michael Kaibel
 */

#ifndef NETWORKIT_BIPARTITE_BIPARTITE_HPP_
#define NETWORKIT_BIPARTITE_BIPARTITE_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

class Bipartite final : public Algorithm {
public:
    /**
     * Creates Bipartite class for graph @G
     *
     * @param G the graph
     */
    Bipartite(const Graph &G);

    /**
     * Calculates if graph @G is bipartite. If yes generates bipartition, if no finds odd cycle
     */
    void run() override;

    /**
     * Check if the graph is bipartite
     *
     * @return true iff the graph is bipartite
     */
    bool isBipartite();

    /**
     * Get a bipartition iff the graph is bipartite
     *
     * @return a bipartition iff the graph is bipartite, throws error otherwise
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

    bool bipartite = false;

    void findOddCircle(std::vector<node> &parent, node v, node w);
};

} // namespace NetworKit

#endif // NETWORKIT_BIPARTITE_BIPARTITE_HPP_

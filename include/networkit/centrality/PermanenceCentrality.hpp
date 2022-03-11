#ifndef NETWORKIT_CENTRALITY_PERMANENCE_CENTRALITY_HPP_
#define NETWORKIT_CENTRALITY_PERMANENCE_CENTRALITY_HPP_

#include <networkit/centrality/Centrality.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 * Permanence centrality measures how well a vertex belongs to its community.
 */
class PermanenceCentrality : public Algorithm {
public:
    /**
     * Constructs the PermanenceCentrality class for the given Graph @a G and Partition @a P.
     *
     * @param G The input graph.
     * @param P Partition for graph G.
     */
    PermanenceCentrality(const Graph &G, const Partition &P);
    void run() override;

    /**
     * Returns the permanence centrality of node @a u.
     *
     * @param u Node in the input graph.
     */
    double getPermanence(node u);

    /**
     * Returns the intra-clustering coefficient of node @a u.
     *
     * @param u Node in the input graph.
     */
    double getIntraClustering(node u);

private:
    const Graph &G;
    const Partition &P;
    std::vector<index> inBegin;
    std::vector<node> inEdges;
    std::vector<bool> marker;
};

} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_PERMANENCE_CENTRALITY_HPP_

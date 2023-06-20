#ifndef NETWORKIT_DISTANCE_DYN_PRUNED_LANDMARK_LABELING_HPP_
#define NETWORKIT_DISTANCE_DYN_PRUNED_LANDMARK_LABELING_HPP_

#include <stdexcept>
#include <vector>

#include <networkit/base/DynAlgorithm.hpp>
#include <networkit/distance/PrunedLandmarkLabeling.hpp>
#include <networkit/dynamics/GraphEvent.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * Dynamic pruned landmark labeling.
 */
class DynPrunedLandmarkLabeling final : public PrunedLandmarkLabeling, public DynAlgorithm {

public:
    /**
     * Dynamic Pruned Landmark Labeling algorithm based on the paper "Fully Dynamic 2-Hop Cover
     * Labeling " from D'Angelo et al., ACM JEA 2019. The algorithm computes distance labels by
     * performing pruned breadth-first searches from each vertex. Distance labels can be updated
     * efficiently after edge insertions.
     * @note this algorithm ignores edge weights and only supports edge insertions.
     *
     * @param G The input graph.
     */
    DynPrunedLandmarkLabeling(const Graph &G) : PrunedLandmarkLabeling(G) {}

    ~DynPrunedLandmarkLabeling() override = default;

    /**
     * Updates the distance labels after an edge insertion on the graph.
     * @note supports only edge insertions.
     *
     * @param e The edge insertion.
     */
    void update(GraphEvent e) override;

    /**
     * Not implemented. The algorithm does not support batch updates.
     * @note This function is not implemented.
     */
    void updateBatch(const std::vector<GraphEvent> &) override {
        throw std::runtime_error("Not implemented.");
    }

private:
    /**
     * Updates the distance labels after an edge insertion.
     * @param u Source node of the edge.
     * @param v Target node of the edge.
     */
    void addEdge(node u, node v);

    void prunedBFS(node k, node startNode, count bfsLevel, bool reverse);

    void sortUpdatedLabels(bool reverse);

    std::vector<node> updatedNodes;
};

} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_DYN_PRUNED_LANDMARK_LABELING_HPP_

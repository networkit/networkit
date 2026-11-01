#ifndef NETWORKIT_DISTANCE_DYN_PRUNED_LANDMARK_LABELING_HPP_
#define NETWORKIT_DISTANCE_DYN_PRUNED_LANDMARK_LABELING_HPP_

#include <span>
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
     * efficiently after edge insertions and edge removals.
     * @note this algorithm ignores edge weights and only supports edge insertions and removals.
     *
     * @param G The input graph.
     */
    DynPrunedLandmarkLabeling(const Graph &G);

    ~DynPrunedLandmarkLabeling() override = default;

    /**
     * Updates the distance labels after an edge insertion or an edge removal on the graph. The
     * graph must already reflect the event.
     * @note supports only edge insertions and edge removals.
     *
     * @param e The edge insertion or edge removal.
     */
    void update(GraphEvent e) override;

    /**
     * Not implemented. The algorithm does not support batch updates.
     * @note This function is not implemented.
     */
    void updateBatch(std::span<const GraphEvent>) override {
        throw std::runtime_error("Not implemented.");
    }

private:
    /**
     * Updates the distance labels after an edge insertion.
     * @param u Source node of the edge.
     * @param v Target node of the edge.
     */
    void addEdge(node u, node v);

    /**
     * Updates the distance labels after an edge removal (DecPLL in the paper): detects the nodes
     * affected by the removal, removes the labels that might be outdated and adds the labels that
     * are needed to restore the cover property.
     * @param x Source node of the edge.
     * @param y Target node of the edge.
     */
    void removeEdge(node x, node y);

    /**
     * Computes the distances from (or to, if @a reverse) @a source.
     */
    std::vector<count> computeDistances(node source, bool reverse) const;

    /**
     * Returns the rank of the smallest hub of the pair (@a u, @a v) at distance @a dist.
     */
    node smallestHub(node u, node v, count dist) const;

    /**
     * Computes the nodes affected by the removal of an edge between @a source and @a target with a
     * BFS from @a source (Algorithm 5 in the paper). The search follows incoming edges if
     * @a reverse is true, and in that case hubs are searched from the visited nodes to @a target.
     * The three distance vectors are relative to the visited nodes.
     *
     * @param source Endpoint of the removed edge the affected nodes are closer to.
     * @param target The other endpoint of the removed edge.
     * @param reverse Whether to follow incoming edges.
     * @param distToTarget Distances from @a target before the removal.
     * @param newDistToTarget Distances from @a target after the removal.
     * @param distToSource Distances from @a source before the removal.
     * @param isAffected Must be all unvisited; set to visited for all affected nodes.
     * @return The affected nodes.
     */
    std::vector<node> detectAffectedNodes(node source, node target, bool reverse,
                                          std::span<const count> distToTarget,
                                          std::span<const count> newDistToTarget,
                                          std::span<const count> distToSource,
                                          std::span<Visit> isAffected);

    /**
     * Removes from the labels of @a nodes all hubs in @a isHub, except for the labels of a node
     * to itself: their distance is always correct (Algorithm 2 in the paper).
     */
    void removeHubs(std::vector<std::vector<Label>> &labels, std::span<const node> nodes,
                    std::span<const Visit> isHub);

    /**
     * Runs a BFS from @a root and adds a label for each visited node in @a isTarget that is not
     * covered at the distance from @a root. The search follows incoming edges if @a reverse is
     * true.
     *
     * @param root Root of the BFS.
     * @param isTarget The nodes to be covered.
     * @param reverse Whether to follow incoming edges.
     * @param ordered If true, nodes with a smaller rank than @a root are skipped (Order-Restore).
     * Otherwise, the hub of a new label is the one with the smaller rank among @a root and the
     * target (Greedy-Restore).
     */
    void restoreCover(node root, std::span<const Visit> isTarget, bool reverse, bool ordered);

    void insertLabel(std::vector<Label> &labelsU, node hub, count dist);

    void prunedBFS(node k, node startNode, count bfsLevel, bool reverse);

    void sortUpdatedLabels(bool reverse);

    std::vector<node> updatedNodes;

    // ranks[u] is the position of u in nodesSortedByDegreeDesc.
    std::vector<node> ranks;
};

} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_DYN_PRUNED_LANDMARK_LABELING_HPP_

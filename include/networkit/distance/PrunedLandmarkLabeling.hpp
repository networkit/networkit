#ifndef NETWORKIT_DISTANCE_PRUNED_LANDMARK_LABELING_HPP_
#define NETWORKIT_DISTANCE_PRUNED_LANDMARK_LABELING_HPP_

#include <utility>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class PrunedLandmarkLabeling : public Algorithm {

public:
    /**
     * Pruned Landmark Labeling algorithm based on the paper "Fast exact shortest-path distance
     * queries on large networks by pruned landmark labeling" from Akiba et al., ACM SIGMOD 2013.
     * The algorithm computes distance labels by performing pruned breadth-first searches from each
     * vertex. Labels are used to quickly retrieve shortest-path distances between node pairs.
     * @note this algorithm only works for unweighted graphs.
     *
     * @param G The input graph.
     */
    PrunedLandmarkLabeling(const Graph &G);

    /**
     * Computes distance labels. Run this function before calling 'query'.
     */
    void run() override;

    /**
     * Returns the shortest-path distance between the two nodes.
     *
     * @param u Source node.
     * @param v Target node.
     *
     * @return The shortest-path distance from @a u to @a v.
     */
    count query(node u, node v) const;

private:
    count queryImpl(node u, node v) const;

    static constexpr count infDist = std::numeric_limits<count>::max();

    struct Label {
        Label() : node_(none), distance_(infDist) {}
        Label(node node_, count distance_) : node_(node_), distance_(distance_) {}
        node node_;
        count distance_;
    };

    const Graph *G;
    std::vector<node> nodesSortedByDegreeDesc;
    std::vector<bool> visited;
    std::vector<std::vector<Label>> labelsOut, labelsIn;
    std::vector<Label> labelsUCopy, labelsVCopy;

    template <bool Reverse = false>
    void prunedBFS(node root, node rankOfRootNode);

    auto getSourceLabelsIterators(node u, bool reverse = false) const {
        if (reverse)
            return std::make_pair(labelsIn[u].begin(), labelsIn[u].end());
        else
            return std::make_pair(labelsOut[u].begin(), labelsOut[u].end());
    }

    auto getTargetLabelsIterators(node u) const {
        return std::make_pair(labelsOut[u].begin(), labelsOut[u].end());
    }
};

} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_PRUNED_LANDMARK_LABELING_HPP_

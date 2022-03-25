/*
 * SPSP.hpp
 *  Created on: 23.10.2020
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_DISTANCE_SPSP_HPP_
#define NETWORKIT_DISTANCE_SPSP_HPP_

#include <unordered_map>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * Class for some-pairs shortest path algorithm.
 */
class SPSP final : public Algorithm {

public:
    /**
     * Creates the SPSP class for @a G.
     *
     * @param G The graph.
     * @param sourcesFirst,sourcesLast Range of the source nodes.
     */
    template <class InputIt>
    explicit SPSP(const Graph &G, InputIt sourcesFirst, InputIt sourcesLast) : G(&G) {
        setSources(sourcesFirst, sourcesLast);
    }

    /**
     * Creates the SPSP class for @a G.
     *
     * @param G The graph.
     * @param sourcesFirst,sourcesLast Range of the source nodes.
     * @param targetsFirst,targetLast Range of the target nodes.
     */
    template <class SourcesInputIt, class TargetsInputIt>
    explicit SPSP(const Graph &G, SourcesInputIt sourcesFirst, SourcesInputIt sourcesLast,
                  TargetsInputIt targetsFirst, TargetsInputIt targetsLast)
        : G(&G) {
        setSources(sourcesFirst, sourcesLast);
        setTargets(targetsFirst, targetsLast);
    }

    /**
     * Sets the source nodes.
     *
     * @param sourcesFirst,sourcesLast Range of the new source nodes.
     */
    template <class InputIt>
    void setSources(InputIt sourcesFirst, InputIt sourcesLast) {
        initSourcesOrTargets(sources, sourcesFirst, sourcesLast, sourceIdx);
    }

    /**
     * Sets the target nodes.
     *
     * @param sourcesFirst,sourcesLast Range of the new source nodes.
     */
    template <class InputIt>
    void setTargets(InputIt targetsFirst, InputIt targetsLast) {
        initSourcesOrTargets(targets, targetsFirst, targetsLast, targetIdx);
    }

    ~SPSP() override = default;

    /**
     * Computes the shortest paths from the source nodes to the target nodes using bidirectional
     * graph explorations. If no targets nodes are provided, the algorithm computes the shortest
     * paths from each source node to all the other nodes. The algorithm is parallel.
     */
    void run() override;

    /**
     * Returns a vector of weighted distances between all the source nodes to all the other nodes.
     *
     * @return The shortest-path distances from the source nodes to any other node in
     * the graph.
     */
    const std::vector<std::vector<edgeweight>> &getDistances() const {
        assureFinished();
        return distances;
    }

    /**
     * Returns the distance from the source @a u to node @a v or infinity if @a
     * u cannot reach @a v.
     *
     * @param u A source node.
     * @param v A node.
     */
    edgeweight getDistance(node u, node v) const {
        assureFinished();
        if (targets.empty())
            return distances[sourceIdx.at(u)][v];
        else
            return distances[sourceIdx.at(u)][targetIdx.at(v)];
    }

    /**
     * Returns the <node, index> map from source nodes to their index in the vector
     * returned by SPSP::getDistances().
     *
     * @return Map from source nodes to their index in the vector returned by SPSP::getDistances().
     */
    const std::unordered_map<node, index> &getSourceIndexMap() const noexcept { return sourceIdx; }

    /**
     * Returns the <node, index> map from target nodes to their index in the vector
     * returned by SPSP::getDistances().
     *
     * @return Map from target nodes to their index in the vector returned by SPSP::getDistances().
     */
    const std::unordered_map<node, index> &getTargetIndexMap() const noexcept { return targetIdx; }

private:
    const Graph *G;
    std::vector<node> sources, targets;
    std::unordered_map<node, index> sourceIdx, targetIdx;
    std::vector<std::vector<edgeweight>> distances;

    template <class InputIt>
    void initSourcesOrTargets(std::vector<node> &vec, InputIt first, InputIt last,
                              std::unordered_map<node, index> &indices) {
        vec.assign(first, last);
        indices.clear();
        for (index i = 0; i < vec.size(); ++i)
            indices.emplace(vec[i], i);
    }

    void runWithTargets();
    void runWithoutTargets();
};

} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_SPSP_HPP_

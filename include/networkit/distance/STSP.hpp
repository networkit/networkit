/*
 * STSP.hpp
 *
 *  Created on: 15.06.2019
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_DISTANCE_STSP_HPP_
#define NETWORKIT_DISTANCE_STSP_HPP_

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

#include <algorithm>
#include <unordered_map>
#include <vector>

namespace NetworKit {

/**
 * @ingroup distance
 * Abstract base class for source-target shortest path algorithms.
 */
class STSP : public Algorithm {

public:
    /**
     * Creates the STSP class for a graph @a G, source node @a source, and
     * target node @a target.
     *
     * @param G The graph.
     * @param source The source node.
     * @param target The target node.
     * @param storePred If true, the algorithm will also store the predecessors
     * and reconstruct a shortest path from @a source and @a target.
     */
    STSP(const Graph &G, node source, node target, bool storePred = true)
        : G(&G), source(source), target(target), storePred(storePred) {
        if (!G.hasNode(source))
            throw std::runtime_error("Error: source node not in the graph!");
        if (!G.hasNode(target))
            throw std::runtime_error("Error: target node not in the graph!");
        if (source == target)
            WARN("Source and target nodes are equal!");
    }

    /**
     * Creates the STSP class for a graph @a G, source node @a source, and
     * multiple target nodes.
     *
     * @param G The graph.
     * @param source The source node.
     * @param targetsFirst,targetsLast Range of target nodes.
     */
    template <class InputIt>
    STSP(const Graph &G, node source, InputIt targetsFirst, InputIt targetsLast)
        : G(&G), source(source), targets(targetsFirst, targetsLast), storePred(false) {
        if (!G.hasNode(source))
            throw std::runtime_error("Error: source node not in the graph!");
        if (!std::all_of(targetsFirst, targetsLast, [&G](node u) { return G.hasNode(u); }))
            throw std::runtime_error("Error: target node not in the graph!");
    }

    /**
     * Returns a shortest path from the source node to the target node (without
     * including them). Note: the shortest path can be constructed only if the
     * algorithm is executed with @a storePred set to true.
     *
     * @return A shortest path from the @a source to @a target.
     */
    virtual std::vector<node> getPath() const {
        checkStorePredecessors();
        assureFinished();
        return path;
    }

    /**
     * Returns the predecessor nodes from the @a target to the source. Note:
     * predecessors are stored only if the algorithm is executed with @a
     * storePred set to true.
     *
     * @return The list of predecessors from @a target to @a source.
     */
    const std::vector<node> &getPredecessors() const {
        checkStorePredecessors();
        assureFinished();
        return pred;
    }

    /**
     * If the target is a single node: returns the distance from the source node to the target
     * node.
     * @return The distance from source to the target node.
     */
    edgeweight getDistance() const {
        assureFinished();
        return distance;
    }

    /**
     * In case of multiple target nodes: returns the distance from the source node to the target
     * nodes.
     * @return Distances from the source to the target nodes.
     */
    const std::vector<edgeweight> &getDistances() const {
        assureFinished();
        return distances;
    }

    /**
     * Sets the source node.
     *
     * @param newSource The new source node.
     */
    void setSource(node newSource) {
        assert(G->hasNode(newSource));
        source = newSource;
    }

    /**
     * Sets the target node.
     *
     * @param newTarget The new target node.
     */
    void setTarget(node newTarget) {
        assert(G->hasNode(newTarget));
        target = newTarget;
        targets.clear();
    }

    /**
     * Sets the target nodes.
     *
     * @param targetsFirst,targetsLast Range of target nodes.
     */
    template <class InputIt>
    void setTargets(InputIt targetsFirst, InputIt targetsLast) {
        assert(std::all_of(targetsFirst, targetsLast, [&](auto u) { return G->hasNode(u); }));
        targets.assign(targetsFirst, targetsLast);
    }

    /**
     * Returns the <node, index> map from target nodes to their index in the vector
     * returned by STSP::getDistances().
     *
     * @return Map from target nodes to their index in the vector returned by STSP::getDistances().
     */
    const std::unordered_map<node, index> &getTargetIndexMap() const {
        assureFinished();
        return targetIdx;
    }

protected:
    const Graph *G;
    node source = none, target = none;
    std::vector<node> targets;
    const bool storePred;
    std::vector<node> pred, path;

    edgeweight distance;
    std::vector<edgeweight> distances;
    std::unordered_map<node, index> targetIdx;

    // Builds a source-target shortest path using the predecessors
    void buildPath();

    void checkStorePredecessors() const {
        if (!storePred)
            WARN("Predecessors not stored.");
    }

    void init() {
        hasRun = false;
        if (storePred)
            pred.resize(G->upperNodeIdBound());
    }
};

inline void STSP::buildPath() {
    path.clear();
    node t = target;
    while (pred[t] != source) {
        path.push_back(pred[t]);
        t = pred[t];
    }
    std::reverse(path.begin(), path.end());
}

} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_STSP_HPP_

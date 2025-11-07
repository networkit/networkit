/*
 * EdgeScore.cpp
 *
 *  Created on: 18.08.2015
 *      Author: Gerd Lindner
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/edgescores/EdgeScore.hpp>

namespace NetworKit {

template <typename T>
EdgeScore<T>::EdgeScore(const Graph &G) : Algorithm(), G(&G), scoreData() {
    if (G.isDirected()) {
        WARN("EdgeScore is not well tested on directed graphs");
    }
}

/** Compute the edge score. */
template <typename T>
void EdgeScore<T>::run() {
    // empty run method for edge scoring methods that do not require preprocessing but calculate
    // score(u,v) on the fly
    hasRun = true;
}

/** Get a vector containing the score for each edge in the graph.
 *
 * @return the edge scores calculated by @link run().
 */
template <typename T>
const std::vector<T> &EdgeScore<T>::scores() const {
    assureFinished();
    return scoreData;
}

/** Get the edge score of the edge with the given edge id.
 */
template <typename T>
T EdgeScore<T>::score(edgeid eid) {
    assureFinished();
    return scoreData[eid];
}

/** Get the edge score of the given edge.
 */
template <typename T>
T EdgeScore<T>::score(node u, node v) {
    return score(G->edgeId(u, v));
}

template <typename T>
Graph EdgeScore<T>::calculate(bool squared, edgeweight offset, edgeweight factor) const {
    assureFinished();

    if (!G->hasEdgeIds()) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }

    // Match existing semantics: weighted, undirected clone of G.
    Graph result(*G, /*weighted=*/true, /*directed=*/false);

    if (squared) {
        G->parallelForEdges([&](node u, node v, edgeid eid) {
            const auto s = scoreData[eid];
            result.setWeight(u, v, offset + factor * s * s);
        });
    } else {
        G->parallelForEdges([&](node u, node v, edgeid eid) {
            result.setWeight(u, v, offset + factor * scoreData[eid]);
        });
    }
    return result;
}

template class EdgeScore<double>;
template class EdgeScore<count>;

} /* namespace NetworKit */

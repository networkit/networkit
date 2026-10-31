/*
 * METISGraphWriterImpl.hpp
 *
 *  Created on: 30.01.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_IO_METIS_GRAPH_WRITER_IMPL_HPP_
#define NETWORKIT_IO_METIS_GRAPH_WRITER_IMPL_HPP_

#include <fstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/io/METISGraphWriter.hpp>

namespace NetworKit {

template <class GraphT>
void GenericMETISGraphWriter<GraphT>::write(const GraphT &G, std::string_view path) {
    this->write(G, G.isWeighted(), path);
}

template <class GraphT>
void GenericMETISGraphWriter<GraphT>::write(const GraphT &G, bool weighted, std::string_view path) {
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

    if (G.isDirected()) {
        throw std::invalid_argument{"METIS does not support directed graphs"};
    }
    std::ofstream file(path.data());
    Aux::enforceOpened(file);

    count n = G.numberOfNodes();
    count m = G.numberOfEdges();

    file << n << " " << m << " " << int{weighted} << '\n';

    // METIS node ids are continuous, so compact them (in node order) if the graph has deleted nodes
    std::vector<count> continuousIds;
    const auto upperNodeIdBound = static_cast<count>(G.upperNodeIdBound());
    if (G.numberOfNodes() != upperNodeIdBound) {
        continuousIds.resize(upperNodeIdBound, none);
        count nextId = 0;
        G.forNodes([&](NodeT u) -> void { continuousIds[u] = nextId++; });
    }
    // METIS-indices are 1-based
    const auto metisId = [&](NodeT v) -> count {
        return (continuousIds.empty() ? static_cast<count>(v) : continuousIds[v]) + 1;
    };

    if (weighted == false) {
        G.forNodes([&](NodeT u) -> void {
            G.forNeighborsOf(u, [&](NodeT v) -> void { file << metisId(v) << " "; });
            file << '\n';
        });
    } else {
        G.forNodes([&](NodeT u) -> void {
            // Shorter callbacks are ambiguous if EdgeWeightT equals NodeT ((v, w) vs. (u, v)) or
            // edgeid ((u, v, w) vs. (u, v, id)). The cast prints 8-bit integral weights as numbers
            // instead of characters.
            G.forNeighborsOf(u, [&](NodeT, NodeT v, EdgeWeightT w, edgeid) -> void {
                file << metisId(v) << " " << static_cast<std::common_type_t<EdgeWeightT, int>>(w)
                     << "\t";
            });
            file << '\n';
        });
    }
}

} /* namespace NetworKit */

#endif // NETWORKIT_IO_METIS_GRAPH_WRITER_IMPL_HPP_

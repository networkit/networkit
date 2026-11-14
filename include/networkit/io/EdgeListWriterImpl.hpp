/*
 * EdgeListWriterImpl.hpp
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#ifndef NETWORKIT_IO_EDGE_LIST_WRITER_IMPL_HPP_
#define NETWORKIT_IO_EDGE_LIST_WRITER_IMPL_HPP_

#include <array>
#include <cassert>
#include <charconv>
#include <fstream>
#include <system_error>
#include <type_traits>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/io/EdgeListWriter.hpp>

namespace NetworKit {

template <class GraphT>
GenericEdgeListWriter<GraphT>::GenericEdgeListWriter(char separator, NodeT firstNode,
                                                     bool bothDirections)
    : separator(separator), firstNode(firstNode), bothDirections(bothDirections) {}

template <class GraphT>
void GenericEdgeListWriter<GraphT>::write(const GraphT &G, std::string_view path) {
    using EdgeWeightT = typename GraphT::EdgeWeightT;

    std::ofstream file(path.data());
    Aux::enforceOpened(file);

    // The sum is computed as count, so that it cannot overflow NodeT. The ids of the file are
    // printed as numbers also for 8-bit node types.
    const auto fileId = [this](NodeT u) -> count {
        return static_cast<count>(u) + static_cast<count>(firstNode);
    };

    // Large enough for the longest text of any arithmetic weight, e.g. -1.7976931348623157e+308.
    std::array<char, 32> buffer;

    // The callbacks take all four arguments, since shorter ones are ambiguous if EdgeWeightT
    // equals NodeT or edgeid.
    auto writeWeightedEdge = [&](NodeT u, NodeT v, EdgeWeightT weight, edgeid) -> void {
        // The cast turns 8-bit integral weights into numbers instead of characters.
        const auto printedWeight = static_cast<std::common_type_t<EdgeWeightT, int>>(weight);
        // std::to_chars writes the shortest text that reads back as the same value, which
        // operator<< does not (it rounds floating-point values to six significant digits).
        [[maybe_unused]] const auto [last, ec] =
            std::to_chars(buffer.data(), buffer.data() + buffer.size(), printedWeight);
        assert(ec == std::errc{});
        file << fileId(u) << separator << fileId(v) << separator;
        file.write(buffer.data(), last - buffer.data());
        file << '\n';
    };

    auto writeUnweightedEdge = [&](NodeT u, NodeT v, EdgeWeightT, edgeid) -> void {
        file << fileId(u) << separator << fileId(v) << '\n';
    };

    if (G.isWeighted()) {
        if (bothDirections) {
            G.forNodes([&](NodeT u) -> void { G.forEdgesOf(u, writeWeightedEdge); });
        } else {
            G.forEdges(writeWeightedEdge);
        }
    } else {
        if (bothDirections) {
            G.forNodes([&](NodeT u) -> void { G.forEdgesOf(u, writeUnweightedEdge); });
        } else {
            G.forEdges(writeUnweightedEdge);
        }
    }

    file.close();
}

} /* namespace NetworKit */

#endif // NETWORKIT_IO_EDGE_LIST_WRITER_IMPL_HPP_

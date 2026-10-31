/*
 * METISGraphReaderImpl.hpp
 *
 *  Created on: 17.01.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_IO_METIS_GRAPH_READER_IMPL_HPP_
#define NETWORKIT_IO_METIS_GRAPH_READER_IMPL_HPP_

#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/StringTools.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/METISParser.hpp>

namespace NetworKit {

template <class GraphT>
GraphT GenericMETISGraphReader<GraphT>::read(std::string_view path) {
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

    METISParser parser(path);

    std::tuple<count, count, index, count> header = parser.getHeader();
    count n = std::get<0>(header);
    count m = std::get<1>(header);
    index fmt = std::get<2>(header);
    count ncon = std::get<3>(header);

    // The node ids 0, ..., n - 1 must be representable by NodeT and differ from NullNodeId<NodeT>,
    // which marks "no node". Since NullNodeId<NodeT> is the largest NodeT, this allows at most
    // max() nodes; the static_assert catches a change of NullNodeId that would invalidate this.
    // The 1-based ids 1, ..., n in the file then fit into NodeT as well.
    static_assert(NullNodeId<NodeT> == std::numeric_limits<NodeT>::max());
    constexpr auto maxNumberOfNodes = static_cast<count>(std::numeric_limits<NodeT>::max());
    if (n > maxNumberOfNodes) {
        throw std::runtime_error("METIS file " + std::string(path) + " has " + std::to_string(n)
                                 + " nodes, but the node type of the graph can hold at most "
                                 + std::to_string(maxNumberOfNodes) + " nodes");
    }

    bool weighted;
    if (fmt % 10 == 1) {
        weighted = true;
        DEBUG("graph has been identified as weighted");
    } else {
        weighted = false;
    }
    count ignoreFirst = 0;
    if (fmt / 10 == 1) {
        DEBUG("first ", ncon, " value(s) will be ignored");
        ignoreFirst = ncon;
    }

    GraphT G(n, weighted);
    std::string graphName =
        Aux::StringTools::split(Aux::StringTools::split(path, '/').back(), '.').front();

    INFO("\n[BEGIN] reading graph G(n=", n, ", m=", m, ") from METIS file: ", graphName);

#ifndef NETWORKIT_RELEASE_LOGGING
    // percentage for progress bar
    double p = 0.0;
#endif // NETWORKIT_RELEASE_LOGGING
    // u is in [0, n) and thus fits into NodeT (checked above)
    count u = 0; // begin with 0
    count edgeCounter = 0;
    count selfLoops = 0;

    // Adds the edge from u to the 1-based METIS node id metisV, returns false if it is skipped.
    const auto addEdge = [&](NodeT metisV, EdgeWeightT weight) -> bool {
        if (metisV == NodeT{0}) {
            ERROR("METIS Node ID should not be 0, edge ignored.");
            return false;
        }
        Aux::Checkers::Enforcer::enforce(metisV > NodeT{0} && static_cast<count>(metisV) <= n);
        const auto v = static_cast<NodeT>(metisV - 1); // METIS-indices are 1-based
        // correct edgeCounter for selfloops
        if (u == static_cast<count>(v)) {
            ++edgeCounter;
            ++selfLoops;
        }
        if (!G.addPartialEdge(unsafe, static_cast<NodeT>(u), v, weight, 0, true)) {
            WARN("Not adding edge ", u, "-", v, " since it is already present.");
            return false;
        }
        return true;
    };

    while (parser.hasNext() && u < n) {
        if (weighted) {
            const std::vector<std::pair<NodeT, EdgeWeightT>> adjacencies =
                parser.getNextWithWeights<NodeT, EdgeWeightT>(ignoreFirst);
            edgeCounter += adjacencies.size();
            DEBUG("node ", u, " has ", adjacencies.size(), " edges");
            for (const auto &[v, weight] : adjacencies) {
                if (addEdge(v, weight)) {
                    // The cast prints 8-bit integral weights as numbers instead of characters.
                    TRACE("(", u, ",", v - 1, ",",
                          static_cast<std::common_type_t<EdgeWeightT, int>>(weight), ")");
                }
            }
        } else {
            const std::vector<NodeT> adjacencies = parser.getNext<NodeT>(ignoreFirst);
            edgeCounter += adjacencies.size();
            for (const NodeT v : adjacencies) {
                addEdge(v, EdgeWeightT{1});
            }
        }
        ++u; // next node
#ifndef NETWORKIT_RELEASE_LOGGING
        if ((u % ((n + 10) / 10)) == 0) {
            p = ((double)(u - 1) / (double)n) * 100;
            DEBUG(p, "% ");
        }
#endif // NETWORKIT_RELEASE_LOGGING
    }
    G.setEdgeCount(unsafe, edgeCounter / 2);
    G.setNumberOfSelfLoops(unsafe, selfLoops);
    if (G.numberOfEdges() != m) {
        ERROR("METIS file ", path,
              " is corrupted: actual number of added edges doesn't match the specifed number of "
              "edges");
    }
    if (edgeCounter != 2 * m) {
        WARN("METIS file is corrupted: not every edge is listed twice");
    }

    INFO("\n[DONE]\n");
    return G;
}

} /* namespace NetworKit */

#endif // NETWORKIT_IO_METIS_GRAPH_READER_IMPL_HPP_

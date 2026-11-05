/*
 * EdgeListReaderImpl.hpp
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#ifndef NETWORKIT_IO_EDGE_LIST_READER_IMPL_HPP_
#define NETWORKIT_IO_EDGE_LIST_READER_IMPL_HPP_

#include <charconv>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <system_error>
#include <type_traits>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/io/EdgeListReader.hpp>
#include <networkit/io/MemoryMappedFile.hpp>

namespace NetworKit {

template <class GraphT>
GenericEdgeListReader<GraphT>::GenericEdgeListReader(char separator, NodeT firstNode,
                                                     std::string_view commentPrefix,
                                                     bool continuous, bool directed)
    : separator(separator), commentPrefix(commentPrefix), firstNode(firstNode),
      continuous(continuous), mapNodeIds(), directed(directed) {
    if (!continuous && firstNode != 0) {
        // firstNode not being 0 in the continuous = false case leads to a segmentation fault
        WARN("firstNode set to 0 since continuous is false");
        this->firstNode = 0;
    }
}

template <class GraphT>
auto GenericEdgeListReader<GraphT>::getNodeMap() const -> const std::map<std::string, NodeT> & {
    if (this->continuous)
        throw std::runtime_error("Input files are assumed to have continuous node ids, therefore "
                                 "no node mapping has been created.");
    return this->mapNodeIds;
}

template <class GraphT>
bool GenericEdgeListReader<GraphT>::atLineEnd(const char *it, const char *end) {
    return it == end || *it == '\r' || *it == '\n';
}

template <class GraphT>
const char *GenericEdgeListReader<GraphT>::skipToLineEnd(const char *it, const char *end) {
    while (!atLineEnd(it, end)) {
        ++it;
    }
    return it;
}

template <class GraphT>
const char *GenericEdgeListReader<GraphT>::skipWhitespaceAndSeparator(const char *it,
                                                                      const char *end) const {
    while (it != end && (*it == ' ' || *it == separator)) {
        ++it;
    }
    return it;
}

template <class GraphT>
const char *GenericEdgeListReader<GraphT>::skipWord(const char *it, const char *end) const {
    while (it != end && *it != ' ' && *it != separator && *it != '\n' && *it != '\r') {
        ++it;
    }
    return it;
}

template <class GraphT>
bool GenericEdgeListReader<GraphT>::isComment(const char *it) const {
    return !commentPrefix.empty() && *it == commentPrefix[0];
}

template <class GraphT>
bool GenericEdgeListReader<GraphT>::isWeighted(const char *it, const char *end) const {
    while (it != end) {
        it = skipWhitespaceAndSeparator(it, end);
        if (!atLineEnd(it, end) && !isComment(it)) {
            // Skip the two node ids.
            for (int i = 0; i < 2; ++i) {
                it = skipWord(it, end);
                it = skipWhitespaceAndSeparator(it, end);
            }
            if (!atLineEnd(it, end)) {
                return true;
            }
        }
        it = skipToLineEnd(it, end);
        if (it != end) {
            ++it;
        }
    }
    return false;
}

template <class GraphT>
GraphT GenericEdgeListReader<GraphT>::read(std::string_view path) {
    using EdgeWeightT = typename GraphT::EdgeWeightT;

    this->mapNodeIds.clear();
    MemoryMappedFile mmfile(path);
    const char *it = mmfile.cbegin();
    const char *end = mmfile.cend();

    DEBUG("separator: ", this->separator);
    DEBUG("first node: ", this->firstNode);

    // This function parses in whole words
    const auto scanWord = [this, &it, &end]() -> std::string {
        const char *start = it;
        it = skipWord(it, end);
        return std::string(start, it);
    };

    const bool weighted = isWeighted(it, end);
    DEBUG("graph has been identified as ", weighted ? "weighted" : "unweighted");

    GraphT graph(0, weighted, directed);

    // Parses the number at it into value and advances it, or throws if it is not a valid T.
    const auto scanNumber = [this, &it, &end](auto &value, const char *what) -> void {
        using T = std::remove_reference_t<decltype(value)>;
        const char *tokenEnd = it;
        while (tokenEnd != end && *tokenEnd != ' ' && *tokenEnd != separator && *tokenEnd != '\n'
               && *tokenEnd != '\r') {
            ++tokenEnd;
        }
        if constexpr (std::is_integral_v<T>) {
            const auto [ptr, ec] = std::from_chars(it, tokenEnd, value);
            // E.g., the fraction of 1.5 is not part of an integer.
            if (ec != std::errc{} || ptr != tokenEnd) {
                throw std::runtime_error(std::string("Scanning ") + what
                                         + " failed. The file may be corrupt.");
            }
            it = ptr;
        } else {
            char *past;
            value = static_cast<T>(std::strtod(it, &past));
            if (past <= it) {
                throw std::runtime_error(std::string("Scanning ") + what
                                         + " failed. The file may be corrupt.");
            }
            it = past;
        }
    };

    const auto scanId = [&, this]() -> NodeT {
        if (continuous) {
            NodeT value;
            scanNumber(value, "node");
            return value;
        } else {
            auto word = scanWord();
            auto found = mapNodeIds.find(word);
            if (found != mapNodeIds.end())
                return found->second;
            if (graph.upperNodeIdBound() == NullNodeId<NodeT>) {
                throw std::runtime_error("The file has more nodes than the node type of the graph "
                                         "can hold");
            }
            auto result = mapNodeIds.insert({word, graph.addNode()});
            if (!result.second)
                throw std::runtime_error("Error in mapping nodes");
            return result.first->second;
        }
    };

    const auto scanWeight = [&]() -> EdgeWeightT {
        EdgeWeightT value;
        scanNumber(value, "weight");
        return value;
    };

    // Map nodes and increase graph size
    const auto mapNode = [&](NodeT in) -> NodeT {
        if (in < firstNode) {
            throw std::runtime_error("Node id " + std::to_string(in)
                                     + " is smaller than the first node "
                                     + std::to_string(firstNode));
        }
        const auto mapped = static_cast<NodeT>(in - firstNode);
        // NullNodeId<NodeT> is the largest NodeT, so that the largest valid node id is smaller.
        if (mapped == NullNodeId<NodeT>) {
            throw std::runtime_error("Node id " + std::to_string(in)
                                     + " does not fit into the node type of the graph");
        }
        if (graph.upperNodeIdBound() <= mapped)
            graph.addNodes(static_cast<count>(mapped - graph.upperNodeIdBound()) + 1);

        return mapped;
    };

    // This function modifies the graph on input.
    const auto handleEdge = [&graph](NodeT source, NodeT target, EdgeWeightT weight) -> void {
        if (!graph.hasEdge(source, target)) {
            graph.addEdge(source, target, weight);
        }
    };

    while (it != end) {
        it = skipWhitespaceAndSeparator(it, end);
        if (!atLineEnd(it, end) && !isComment(it)) {
            auto sourceId = scanId();
            if (it == end)
                throw std::runtime_error("Unexpected end of file");
            if (!(*it == ' ' || *it == separator))
                throw std::runtime_error(
                    "Error in parsing file - pointer is whitespace or separator");
            it = skipWhitespaceAndSeparator(it, end);

            auto targetId = scanId();
            it = skipWhitespaceAndSeparator(it, end);

            // Edges without a weight in a weighted graph get the default weight.
            EdgeWeightT edgeWeight{1};
            if (weighted && !atLineEnd(it, end)) {
                edgeWeight = scanWeight();
            }
            handleEdge(mapNode(sourceId), mapNode(targetId), edgeWeight);
        }
        it = skipToLineEnd(it, end);
        if (it != end) {
            ++it;
        }
    }

    graph.shrinkToFit();
    return graph;
}

} /* namespace NetworKit */

#endif // NETWORKIT_IO_EDGE_LIST_READER_IMPL_HPP_

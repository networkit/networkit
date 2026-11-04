/*
 * EdgeListReader.hpp
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#ifndef NETWORKIT_IO_EDGE_LIST_READER_HPP_
#define NETWORKIT_IO_EDGE_LIST_READER_HPP_

#include <map>
#include <string>

#include <networkit/io/GraphReader.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * A reader for various edge list formats, in which each line contains an edge as
 * two node ids.
 */
template <class GraphT>
class GenericEdgeListReader final : public GraphReaderBase<GraphT> {
    using NodeT = typename GraphT::NodeT;

public:
    GenericEdgeListReader() = default; // nullary constructor for Python shell

    /**
     * @param[in]  separator  character used to separate nodes in an edge line
     * @param[in]  firstNode  index of the first node in the file
     * @param[in]  commentPrefix  prefix of comment lines
     * @param[in]  continuous  boolean to specify if node ids are continuous
     * @param[in]  directed  read graph as directed
     *
     * The graph is weighted iff at least one edge line has a third entry, the edge weight; edges
     * without a weight then get the default weight 1. Reading fails if a node id or the number of
     * nodes does not fit into GraphT::NodeT, or if GraphT::EdgeWeightT is integral and a weight is
     * not an integer representable by it.
     */
    GenericEdgeListReader(char separator, NodeT firstNode, std::string_view commentPrefix = "#",
                          bool continuous = true, bool directed = false);

    /**
     * Given the path of an input file, read the graph contained.
     *
     * @param[in]  path  input file path
     */
    GraphT read(std::string_view path) override;

    /**
     * Return the node map, in case node ids are not continuous
     */
    const std::map<std::string, NodeT> &getNodeMap() const;

private:
    // Whether it points to a line ending or to the end of the file.
    static bool atLineEnd(const char *it, const char *end);
    // Returns the position of the next line ending at or after it, or end.
    static const char *skipToLineEnd(const char *it, const char *end);
    // Returns the position of the first character at or after it that is no space or separator.
    const char *skipWhitespaceAndSeparator(const char *it, const char *end) const;
    // Returns the end of the word at it, i.e., the next space, separator, or line ending.
    const char *skipWord(const char *it, const char *end) const;
    // Whether the line at it is a comment. Only call it with it != end.
    bool isComment(const char *it) const;
    // Sweeps over [begin, end): the graph is weighted iff at least one edge line has a third entry.
    bool isWeighted(const char *begin, const char *end) const;

    char separator; //!< character separating nodes in an edge line
    std::string commentPrefix;
    NodeT firstNode;
    bool continuous;
    std::map<std::string, NodeT> mapNodeIds;
    bool directed;
};

using EdgeListReader = GenericEdgeListReader<Graph>;

} /* namespace NetworKit */

#include <networkit/io/EdgeListReaderImpl.hpp>

#endif // NETWORKIT_IO_EDGE_LIST_READER_HPP_

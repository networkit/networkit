/*
 * GraphEvent.hpp
 *
 *  Created on: 02.04.2013
 *      Author: cls
 */

#ifndef NETWORKIT_DYNAMICS_GRAPH_EVENT_HPP_
#define NETWORKIT_DYNAMICS_GRAPH_EVENT_HPP_

#include <networkit/graph/Graph.hpp>

#include <tlx/define/deprecated.hpp>

namespace NetworKit {

/**
 * @ingroup dynamics
 */
class GraphEvent final {

public:
    enum Type {
        NODE_ADDITION,
        NODE_REMOVAL,
        NODE_RESTORATION,
        EDGE_ADDITION,
        EDGE_REMOVAL,
        EDGE_WEIGHT_UPDATE,
        EDGE_WEIGHT_INCREMENT,
        TIME_STEP
    };

    Type type;    //!< type of graph event
    node u;       //!< first node parameter
    node v;       //!< second node parameter
    edgeweight w; //!< edge weight parameter

    GraphEvent() = default;

    GraphEvent(Type type, node u = none, node v = none, edgeweight w = defaultEdgeWeight);

    static bool TLX_DEPRECATED(compare(GraphEvent a, GraphEvent b));
    static bool TLX_DEPRECATED(equal(GraphEvent a, GraphEvent b));

    bool operator==(const GraphEvent &rhs) const noexcept;
    bool operator!=(const GraphEvent &rhs) const noexcept;
    friend bool operator<(const GraphEvent &a, const GraphEvent &b) noexcept;
    friend bool operator>(const GraphEvent &a, const GraphEvent &b) noexcept;
    bool operator<=(const GraphEvent &rhs) const noexcept;
    bool operator>=(const GraphEvent &rhs) const noexcept;

    /**
     * Return string representation.
     */
    std::string toString() const;
};

} /* namespace NetworKit */
#endif // NETWORKIT_DYNAMICS_GRAPH_EVENT_HPP_

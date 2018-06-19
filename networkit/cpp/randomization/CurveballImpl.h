/*
 * CurveballImpl.h
 *
 * Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */
#ifndef RANDOMIZATION_CURVEBALL_IMPL_H
#define RANDOMIZATION_CURVEBALL_IMPL_H

#include <cassert>
#include <utility>

#include "../Globals.h"
#include "../graph/Graph.h"

namespace NetworKit {
namespace CurveballDetails {

// Global Definitions
using TradeDescriptor = std::pair<node, node>;
using tradeid_t = node;
using count = node;

using trade_vector = std::vector<TradeDescriptor>;
using nodepair_vector = std::vector <std::pair<node, node> >;

constexpr node INVALID_NODE = std::numeric_limits<node>::max();
constexpr count LISTROW_END = std::numeric_limits<count>::max();
constexpr tradeid_t TRADELIST_END = std::numeric_limits<tradeid_t>::max();

struct edge_t : public std::pair<node, node> {
    edge_t() : std::pair<node, node>() {}

    edge_t(const std::pair <node, node> &edge) : std::pair<node, node>(
        edge) {}

    edge_t(const node &v1, const node &v2) : std::pair<node, node>(v1, v2) {}

    static edge_t invalid() {
        return edge_t(INVALID_NODE, INVALID_NODE);
    }

    void normalize() {
        if (first > second)
            std::swap(first, second);
    }

    bool is_invalid() const {
        return (*this == invalid());
    }
};


class CurveballAdjacencyList {
public:
    using degree_vector = std::vector<count>;
    using neighbour_vector = std::vector<node>;
    using pos_vector = std::vector<edgeid>;
    using pos_it = pos_vector::iterator;
    using neighbour_it = neighbour_vector::iterator;
    using cneighbour_it = neighbour_vector::const_iterator;
    using nodepair_vector = std::vector< std::pair<node, node> >;

protected:

    neighbour_vector _neighbours;
    degree_vector _offsets;
    pos_vector _begin;
    edgeid _degree_count;

public:
    CurveballAdjacencyList() = default;

    // Receives the degree_vector to initialize
    // As trades permute neighbours the degrees don't change
    CurveballAdjacencyList(const degree_vector& degrees, const edgeid degree_count);

    void initialize(const degree_vector& degrees, const edgeid degree_count);

    void restructure();

    // No Copy Constructor
    CurveballAdjacencyList(const CurveballAdjacencyList &) = delete;

    neighbour_it begin(const node node_id);

    neighbour_it end(const node node_id);

    cneighbour_it cbegin(const node node_id) const;

    cneighbour_it cend(const node node_id) const;

    nodepair_vector getEdges() const;

    void insertNeighbour(const node node_id, const node neighbour) {
        auto pos = begin(node_id) + _offsets[node_id];

        assert(*pos != LISTROW_END);

        *pos = neighbour;

        _offsets[node_id]++;
    }

    node numberOfNodes() const {
        return static_cast<node>(_offsets.size());
    }

    node numberOfEdges() const {
        return static_cast<edgeid>(_degree_count);
    }

    void resetRow(const node node_id) {
        assert(node_id >= 0);
        assert(node_id < static_cast<node>(_offsets.size()));

        _offsets[node_id] = 0;

        return;
    }

    count degreeAt(node node_id) const {
        assert(node_id < static_cast<node>(_offsets.size()));
        assert(node_id >= 0);

        return _begin[node_id + 1] - _begin[node_id] - 1;
    }
};


class CurveballMaterialization {

protected:
    const CurveballAdjacencyList& _adj_list;

public:
    CurveballMaterialization(const CurveballAdjacencyList& adj_list);

    Graph toGraph(bool parallel);

protected:
    void toGraphParallel(Graph &G);
    void toGraphSequential(Graph &G);
};

class TradeList {
public:
    using edge_vector = std::vector<edge_t>;
    using offset_vector = std::vector<tradeid_t>;
    using tradeid_vector = std::vector<tradeid_t>;
    using trade = TradeDescriptor;
    using trade_vector = std::vector<trade>;
    using tradeid_it = tradeid_vector::const_iterator;

protected:
    tradeid_vector _trade_list;
    offset_vector _offsets;
    const node _num_nodes;

public:
    TradeList(const node num_nodes);

    // Receives the edge_vector to initialize
    TradeList(const trade_vector& trades, const node num_nodes);

    // Initialize method
    void initialize(const trade_vector& trades);

    // No Copy Constructor
    TradeList(const TradeList&) = delete;

    tradeid_it getTrades(const node nodeid) const {
        assert(nodeid >= 0);
        assert(nodeid < _num_nodes);

        return _trade_list.begin() + _offsets[nodeid];
    }

    void incrementOffset(const node nodeid) {
        assert(nodeid >= 0);
        assert(nodeid < _num_nodes);
        assert(1 <= _offsets[nodeid + 1] - _offsets[nodeid]);

        _offsets[nodeid]++;
    }

    node numberOfNodes() const {
        return _num_nodes;
    }
};

class CurveballIM  {
public:
    CurveballIM(const Graph& G);

    void run(const trade_vector& trades);

    count getNumberOfAffectedEdges() const {
        assert(_hasRun);
        return _aff_edges;
    }

    Graph getGraph(bool parallel) const;

    nodepair_vector getEdges() const;


protected:
    const Graph& _G;
    const node _num_nodes;

    bool _hasRun;
    CurveballAdjacencyList _adj_list;
    TradeList _trade_list;
    count _max_degree;
    edgeid _aff_edges; // affected half-edges

    void load_from_graph(const trade_vector& trades);

    void restructure_graph(const trade_vector& trades);

    inline void update(const node a, const node b) {
        const tradeid_t ta = *(_trade_list.getTrades(a));
        const tradeid_t tb = *(_trade_list.getTrades(b));
        if (ta < tb) {
            _adj_list.insertNeighbour(a, b);
            return;
        }

        if (ta > tb) {
            _adj_list.insertNeighbour(b, a);
            return;
        }
        // ta == tb
        {
            _adj_list.insertNeighbour(a, b);
        }
    }

};


} // namespace CurveballDetails
} // namespace NetworKit

#endif // ! RANDOMIZATION_CURVEBALL_IMPL_H
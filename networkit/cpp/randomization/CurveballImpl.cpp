/*
 * CurveballImpl.cpp
 *
 * Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */

#include "CurveballImpl.h"

#include <cassert>
#include <algorithm>
#include <vector>

#include "../auxiliary/Timer.h"
#include "../auxiliary/SignalHandling.h"
#include "../auxiliary/RandomBipartitionShuffle.h"


namespace NetworKit {
namespace CurveballDetails {

using degree_vector = std::vector<count>;
using trade_vector = std::vector<TradeDescriptor>;
using neighbour_vector = std::vector<node>;
using node_vector = std::vector<node>;
using nodepair_vector = std::vector< std::pair<node, node> >;


using neighbour_vector = std::vector<node>;
using degree_vector = std::vector<count>;
using degree_it = std::vector<count>::const_iterator;
using pos_vector = std::vector<edgeid>;
using neighbour_it = neighbour_vector::iterator;
using cneighbour_it = neighbour_vector::const_iterator;
using nodepair_vector = std::vector< std::pair<node, node> >;

///////////////////////////////////////////////////////////////////////////////

// public static constexpr count LISTROW_END = std::numeric_limits<count>::max();
/**
 * @brief Initialize method (when constructor can't be used)
 *
 */
void CurveballAdjacencyList::initialize(const degree_vector& degrees,
                                        const edgeid degree_count) {
    _neighbours.resize(degree_count + degrees.size() + 1);
    _offsets.resize(degrees.size());
    _begin.resize(degrees.size() + 1);
    _degree_count = degree_count;

    count sum = 0;
    node node_id = 0;
    for (const count node_degree : degrees) {
        _begin[node_id] = sum;

        assert(node_degree > 0);

        sum += node_degree;
        _neighbours[sum] = LISTROW_END;

        // shift after Sentinel
        sum += 1;
        node_id++;
    }
    _neighbours[sum] = LISTROW_END;
    _begin[degrees.size()] = sum;

    assert(sum == degree_count + degrees.size());
    assert(node_id == degrees.size());

    return;
}

void CurveballAdjacencyList::restructure() {
    std::fill(_offsets.begin(), _offsets.end(), 0);
    return;
}

/**
 * @brief Constructor
 * @param degree_vector Pointer to a vector with node degrees
 * @param degree_count Sum of all degrees in degree_vector
 *
 * We add to each adjacency list entry a delimiter to mark the end
 */
CurveballAdjacencyList::CurveballAdjacencyList(const degree_vector& degrees,
                                               const edgeid degree_count)
    : _neighbours(degree_count + degrees.size() + 1)
    , _offsets(degrees.size())
    , _begin(degrees.size() + 1)
    , _degree_count(degree_count)
{
    count sum = 0;
    node node_id = 0;
    for (const count node_degree : degrees) {
        _begin[node_id] = sum;

        // no isolated nodes allowed
        assert(node_degree > 0);

        sum += node_degree;
        _neighbours[sum] = LISTROW_END;

        sum += 1;
        node_id++;
    }
    _neighbours[sum] = LISTROW_END;
    _begin[degrees.size()] = sum;

    assert(sum == degree_count + degrees.size());
    assert(node_id == degrees.size());
}

neighbour_it CurveballAdjacencyList::begin(const node node_id) {
    return _neighbours.begin() + _begin[node_id];
}

neighbour_it CurveballAdjacencyList::end(const node node_id) {
    return _neighbours.begin() + _begin[node_id] + _offsets[node_id];
}

cneighbour_it CurveballAdjacencyList::cbegin(const node node_id) const {
    return _neighbours.cbegin() + _begin[node_id];
}

cneighbour_it CurveballAdjacencyList::cend(const node node_id) const {
    return _neighbours.cbegin() + _begin[node_id] + _offsets[node_id];
}

nodepair_vector CurveballAdjacencyList::getEdges() const {
    nodepair_vector edges;
    edges.reserve(_degree_count);
    for (node nodeid = 0; nodeid < static_cast<node>(_offsets.size()); nodeid++) {
        for (auto it = cbegin(nodeid); it != cend(nodeid); it++) {
            edges.push_back(std::make_pair(nodeid, *it));
        }
    }

    return edges;
}

///////////////////////////////////////////////////////////////////////////////

CurveballMaterialization::CurveballMaterialization(const CurveballAdjacencyList& adj_list)
    : _adj_list(adj_list)
{ }

NetworKit::Graph CurveballMaterialization::toGraph(bool parallel) {
    Graph G(_adj_list.numberOfNodes(), false, false);

    if (parallel)
        toGraphParallel(G);
    else
        toGraphSequential(G);

    return G;
}

void CurveballMaterialization::toGraphParallel(Graph &G) {
    // 1) insertion of first half is threadsafe
    // 2) insertion of second half is not threadsafe, for now done sequentially

    const node numNodes = _adj_list.numberOfNodes();

    std::vector<NetworKit::count> missingEdgesCounts(numNodes, 0);

    std::vector < std::vector<edgeweight> > new_edgeWeights(numNodes);
    std::vector < std::vector<node> > new_outEdges(numNodes);

    // Add first half of edges and count missing edges for each node
    G.parallelForNodes([&](node nodeid) {
        const count degree =
            static_cast<count>(_adj_list.cend(nodeid) - _adj_list.cbegin(nodeid));
        G.outDeg[nodeid] = degree;
        missingEdgesCounts[nodeid] = _adj_list.degreeAt(nodeid) - degree;
        new_outEdges[nodeid].reserve(degree);
        new_edgeWeights[nodeid].resize(degree, 1);
        for (auto it = _adj_list.cbegin(nodeid); it != _adj_list.cend(nodeid); it++) {
            new_outEdges[nodeid].push_back(*it);
        }
    });

    G.outEdges.swap(new_outEdges);

    // Reserve the space
    G.parallelForNodes([&](node v) {
        G.outEdges[v].reserve(G.outDeg[v] + missingEdgesCounts[v]);
    });

    // Second half of the edges
    G.forNodes([&](node v) {
        for (count neighbor_id = 0; neighbor_id < G.outDeg[v]; neighbor_id++) {
            const node u = G.outEdges[v][neighbor_id];
            G.outEdges[u].push_back(v);
        }
    });

    //TODO: is the networkit adjacency list even sorted for the neighbours? if not omit this
    // Sort neighbours
    G.parallelForNodes([&](node v) {
        std::sort(G.outEdges[v].begin(), G.outEdges[v].end());
    });

    // Set node degrees
    G.parallelForNodes([&](node v) {
        G.outDeg[v] = _adj_list.degreeAt(v);
    });

    // Set number of self-loops
    G.storedNumberOfSelfLoops = 0;

    // Set numberOfEdges
    G.m = _adj_list.numberOfEdges() / 2;

    // Shrink to fit
    G.shrinkToFit();
}

void CurveballMaterialization::toGraphSequential(Graph &G) {
    const node numNodes = _adj_list.numberOfNodes();

    // Analogue to "toGraphSequential" of GraphBuilder
    std::vector<NetworKit::count> missingEdgesCounts;
    missingEdgesCounts.reserve(numNodes);

    std::vector < std::vector<edgeweight> > new_edgeWeights(numNodes);
    std::vector < std::vector<node> > new_outEdges(numNodes);

    // Add first half of edges and count missing edges for each node
    G.forNodes([&](node nodeid) {
        const count degree =
            static_cast<count>(_adj_list.cend(nodeid) - _adj_list.cbegin(nodeid));
        G.outDeg[nodeid] = degree;
        missingEdgesCounts.push_back(_adj_list.degreeAt(nodeid) - degree);
        new_outEdges[nodeid].reserve(degree);
        new_edgeWeights[nodeid].resize(degree, 1);
        for (auto it = _adj_list.cbegin(nodeid); it != _adj_list.cend(nodeid); it++) {
            new_outEdges[nodeid].push_back(*it);
        }
    });

    G.outEdges.swap(new_outEdges);

    // Reserve the space
    G.forNodes([&](node v) {
        G.outEdges[v].reserve(G.outDeg[v] + missingEdgesCounts[v]);
    });

    // Second half of the edges
    G.forNodes([&](node v) {
        for (count neighbor_id = 0; neighbor_id < G.outDeg[v]; neighbor_id++) {
            const node u = G.outEdges[v][neighbor_id];
            G.outEdges[u].push_back(v);
        }
    });

    //TODO: is the networkit adjacency list even sorted for the neighbours? if not omit this
    // Sort neighbours
    G.forNodes([&](node v) {
        std::sort(G.outEdges[v].begin(), G.outEdges[v].end());
    });

    // Set node degrees
    G.forNodes([&](node v) {
        G.outDeg[v] = _adj_list.degreeAt(v);
    });

    // Set number of self-loops
    G.storedNumberOfSelfLoops = 0;

    // Set numberOfEdges
    G.m = _adj_list.numberOfEdges() / 2;

    // Shrink to fit
    G.shrinkToFit();
}

///////////////////////////////////////////////////////////////////////////////
TradeList::TradeList(const node num_nodes)
    : _num_nodes(num_nodes)
{ }

void TradeList::initialize(const trade_vector& trades) {
    _trade_list.clear();
    _trade_list.resize(2 * trades.size() + _num_nodes);
    _offsets.clear();
    _offsets.resize(_num_nodes + 1);

    assert(_num_nodes > 0);
    assert(trades.size() > 0);

    std::vector<tradeid_t> trade_count(_num_nodes);

    // Push occurrences
    for (const auto trade : trades) {
        assert(trade.first >= 0);
        assert(trade.first < _num_nodes);
        assert(trade.second >= 0);
        assert(trade.second < _num_nodes);

        trade_count[trade.first]++;
        trade_count[trade.second]++;
    }

    // add missing +1 for sentinel
    trade_count[0]++;
    std::partial_sum(trade_count.cbegin(), trade_count.cend(), _offsets.begin() + 1, [&](const tradeid_t a, const tradeid_t b){
        return a + b + 1;
    });
    // add dummy
    _offsets[_num_nodes] = 2 * trades.size() + _num_nodes - 1;

    // set sentinels
    for (node node = 1; node < _num_nodes; node++) {
        _trade_list[_offsets[node] - 1] = TRADELIST_END;
    }
    // set last entry as sentinel
    _trade_list.back() = TRADELIST_END;

    std::fill(trade_count.begin(), trade_count.end(), 0);
    {
        tradeid_t trade_id = 0;
        for (const auto trade : trades) {
            auto updateNode = [&] (const node u) {
                const node pos = _offsets[u] + trade_count[u];
                _trade_list[pos] = trade_id;
                trade_count[u]++;
            };

            updateNode(trade.first);
            updateNode(trade.second);
            trade_id++;
        }
    }
}

TradeList::TradeList(const trade_vector& trades, const node num_nodes)
    : _trade_list(2 * trades.size() + num_nodes)
    , _offsets(num_nodes + 1)
    , _num_nodes(num_nodes)
{
    // Manuel: see above

    assert(num_nodes > 0);
    assert(trades.size() > 0);

    std::vector<tradeid_t> trade_count(num_nodes);

    // Push occurences
    for (const auto trade : trades) {
        assert(trade.first >= 0);
        assert(trade.first < num_nodes);
        assert(trade.second >= 0);
        assert(trade.second < num_nodes);

        trade_count[trade.first]++;
        trade_count[trade.second]++;
    }

    // add missing +1 for sentinel
    trade_count[0]++;
    std::partial_sum(trade_count.cbegin(), trade_count.cend(), _offsets.begin() + 1, [&](const tradeid_t a, const tradeid_t b){
        return a + b + 1;
    });
    // add dummy
    _offsets[num_nodes] = 2 * trades.size() + num_nodes - 1;

    // set sentinels
    for (node node = 1; node < _num_nodes; node++) {
        _trade_list[_offsets[node] - 1] = TRADELIST_END;
    }
    // set last entry as sentinel
    _trade_list.back() = TRADELIST_END;

    std::fill(trade_count.begin(), trade_count.end(), 0);
    {
        tradeid_t trade_id = 0;
        for (const auto trade : trades) {
            auto updateNode = [&] (const node u) {
                const node pos = _offsets[u] + trade_count[u];
                _trade_list[pos] = trade_id;
                trade_count[u]++;
            };

            updateNode(trade.first);
            updateNode(trade.second);
            trade_id++;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////




CurveballIM::CurveballIM(const NetworKit::Graph& G)
    : _G(G)
    , _num_nodes(G.numberOfNodes())  
    , _trade_list(G.numberOfNodes())
    , _aff_edges(0)
{
    _hasRun = false;
    assert(G.checkConsistency());
    assert(G.numberOfSelfLoops() == 0);
    assert(_num_nodes > 0);
}

void CurveballIM::load_from_graph(const trade_vector& trades) {
    // Compute degree sequence
    degree_vector degrees;
    degrees.reserve(_num_nodes);
    edgeid degree_sum = 0;
    _G.forNodes([&](node v) {
        degrees.push_back(_G.degree(v));
        degree_sum += _G.degree(v);
    });

    _max_degree = *(std::max_element(degrees.cbegin(), degrees.cend()));

    _adj_list.initialize(degrees, degree_sum);
    _trade_list.initialize(trades);

    // Insert to adjacency list, directed according trades
    _G.forEdges([&](node u, node v) {
        update(u, v);
    });
    return;
}

void CurveballIM::restructure_graph(const trade_vector& trades) {
    nodepair_vector edges =_adj_list.getEdges();

    _adj_list.restructure();
    _trade_list.initialize(trades);

    for (const auto edge : edges) {
        update(edge.first, edge.second);
    }

    return;
}

void CurveballIM::run(const trade_vector& trades) {
    if (!_hasRun)
        load_from_graph(trades);
    else
        restructure_graph(trades);

    NetworKit::count trade_count = 0;
    neighbour_vector common_neighbours;
    neighbour_vector disjoint_neighbours;

    common_neighbours.reserve(_max_degree);
    disjoint_neighbours.reserve(_max_degree);

    Aux::SignalHandler handler;

    auto& urng = Aux::Random::getURNG();

    for (const auto& trade : trades) {
        handler.assureRunning();

        // Trade partners u and v
        const node u = trade.first;
        const node v = trade.second;

        _aff_edges += _adj_list.degreeAt(u);
        _aff_edges += _adj_list.degreeAt(v);

        // Shift the _trade_list offset for these two, currently was set to trade_count
        _trade_list.incrementOffset(u);
        _trade_list.incrementOffset(v);

        // Retrieve respective neighbours
        // we return whether u has v in his neighbors or vice-versa
        auto organize_neighbors = [&](node node_x, node node_y) {
            auto pos = std::find(_adj_list.begin(node_x), _adj_list.end(node_x), node_y);
            if (pos == _adj_list.cend(node_x)) {
                // element not found, sort anyway
                std::sort(_adj_list.begin(node_x), _adj_list.end(node_x));

                return false;
            } else {
                // overwrite node_y's position with END
                *pos = LISTROW_END;

                // sort, such that node_y's position is at end - 1
                std::sort(_adj_list.begin(node_x), _adj_list.end(node_x));

                // overwrite with node_y again
                *(_adj_list.end(node_x) - 1) = node_y;

                return true;
            }
        };

        const bool u_share = organize_neighbors(u, v);
        const bool v_share = organize_neighbors(v, u);
        auto u_end = (u_share ? _adj_list.cend(u) - 1 : _adj_list.cend(u));
        auto v_end = (v_share ? _adj_list.cend(v) - 1 : _adj_list.cend(v));

        const bool shared = u_share || v_share;

        // both can't have each other, only inserted in one
        assert((!u_share && !v_share) || (u_share != v_share));

        // No need to keep track of direct positions
        // Get common and disjoint neighbors
        // Here sort and parallel scan
        common_neighbours.clear();
        disjoint_neighbours.clear();
        auto u_nit = _adj_list.cbegin(u);
        auto v_nit = _adj_list.cbegin(v);
        while ((u_nit != u_end) && (v_nit != v_end)) {
            assert(*u_nit != v);
            assert(*v_nit != u);
            if (*u_nit > *v_nit) {
                disjoint_neighbours.push_back(*v_nit);
                v_nit++;
                continue;
            }
            if (*u_nit < *v_nit) {
                disjoint_neighbours.push_back(*u_nit);
                u_nit++;
                continue;
            }
            // *u_nit == *v_nit
            {
                common_neighbours.push_back(*u_nit);
                u_nit++;
                v_nit++;
            }
        }
        if (u_nit == u_end)
            disjoint_neighbours.insert(disjoint_neighbours.end(), v_nit, v_end);
        else
            disjoint_neighbours.insert(disjoint_neighbours.end(), u_nit, u_end);

        const count u_setsize = static_cast<count>(u_end - _adj_list.cbegin(u) - common_neighbours.size());
        const count v_setsize = static_cast<count>(v_end - _adj_list.cbegin(v) - common_neighbours.size());
        // v_setsize not necessarily needed

        // Reset fst/snd row
        _adj_list.resetRow(u);
        _adj_list.resetRow(v);

        Aux::random_bipartition_shuffle(disjoint_neighbours.begin(), disjoint_neighbours.end(),
                         u_setsize, urng);

        // Assign first u_setsize to u and last v_setsize to v
        // if not existent then max value, and below compare goes in favor of partner, if partner
        // has no more neighbours as well then their values are equal (max and equal)
        // and tiebreaking is applied
        for (count counter = 0; counter < u_setsize; counter++) {
            const node swapped = disjoint_neighbours[counter];
            update(u, swapped);
        }
        for (count counter = u_setsize; counter < u_setsize + v_setsize; counter++) {
            const node swapped = disjoint_neighbours[counter];
            update(v, swapped);
        }
        // Distribute common edges
        for (const auto common : common_neighbours) {
            update(u, common);
            update(v, common);
        }
        // Do not forget edge between u and v
        if (shared)
            update(u, v);

        trade_count++;
    }

    _hasRun = true;

    return;
}

NetworKit::Graph CurveballIM::getGraph(bool parallel) const {
    CurveballMaterialization gb(_adj_list);

    return gb.toGraph(parallel);
}


} // namespace CurveballDetails
} // namespace NetworKit

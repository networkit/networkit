#include <networkit/graph/DHBGraph.hpp>

#include <cmath>
#include <map>
#include <random>
#include <sstream>
#include <thread>

namespace NetworKit {

/** CONSTRUCTORS **/

DHBGraph::DHBGraph(count n, bool weighted, bool directed, bool edgesIndexed)
    : weighted(weighted), directed(directed), m(0), storedNumberOfSelfLoops(0), omega(0),
      edgesIndexed(edgesIndexed), // edges are not indexed by default
      nodeAttributeMap(this), edgeAttributeMap(this) {

    m_dhb_graph = dhb::Matrix<EdgeData>(n);
}

DHBGraph::DHBGraph(std::initializer_list<WeightedEdge> edges) : DHBGraph(0, true) {
    using namespace std;

    m_dhb_graph = dhb::Matrix<EdgeData>();
    for (const auto &weighted_edge : edges) {
        auto maxId = std::max(weighted_edge.u, weighted_edge.v);
        if (maxId >= m_dhb_graph.vertices_count()) {
            m_dhb_graph.resize(maxId + 1);
        }
        addEdge(weighted_edge.u, weighted_edge.v, weighted_edge.weight);
    }
}

/** EDGE IDS **/

void DHBGraph::indexEdges(bool force) {
    if (edgesIndexed && !force) {
        return;
    }

    omega = 0; // reset edge ids (for re-indexing)

    if (directed) {
        forEdges([&](node u, node v) {
            m_dhb_graph.neighbors(u).iterator_to(v)->data().id = omega;
            ++omega;
        });
    } else {
        forEdges([&](node u, node v) {
            if (u == v) { // first time encounter self-loop
                m_dhb_graph.neighbors(u).iterator_to(v)->data().id = omega;
                ++omega;
            } else {
                m_dhb_graph.neighbors(u).iterator_to(v)->data().id = omega;
                m_dhb_graph.neighbors(v).iterator_to(u)->data().id = omega;
                ++omega;
            }
        });
    }
    edgesIndexed = true;
}

edgeid DHBGraph::edgeId(node u, node v) const {
    if (!edgesIndexed) {
        throw std::runtime_error(
            "Edges have not been indexed! To fix this, simply call indexEdges first.");
    }
    return m_dhb_graph.neighbors(u).iterator_to(v)->data().id;
}

/** GRAPH INFORMATION **/

bool DHBGraph::isIsolated(node v) const {
    if (directed) {
        bool const has_no_outgoing_edges = (m_dhb_graph.degree(v) == 0);
        bool has_no_incoming_edges = true;
        forEdges([&](node, node to) {
            if (v == to) {
                has_no_incoming_edges = false;
            }
        });
        return has_no_outgoing_edges && has_no_incoming_edges;
    }
    return (m_dhb_graph.degree(v) == 0);
}

edgeweight DHBGraph::computeWeightedDegree(node u, bool inDegree, bool countSelfLoopsTwice) const {
    if (weighted) {
        edgeweight sum = 0.0;
        auto sumWeights = [&](node v, edgeweight w) {
            sum += (countSelfLoopsTwice && u == v) ? 2. * w : w;
        };
        if (inDegree) {
            forInNeighborsOf(u, sumWeights);
        } else {
            forNeighborsOf(u, sumWeights);
        }
        return sum;
    }

    count sum = inDegree ? degreeIn(u) : degreeOut(u);
    auto countSelfLoops = [&](node v) { sum += (u == v); };

    if (countSelfLoopsTwice && numberOfSelfLoops()) {
        if (inDegree) {
            forInNeighborsOf(u, countSelfLoops);
        } else {
            forNeighborsOf(u, countSelfLoops);
        }
    }

    return static_cast<edgeweight>(sum);
}

/** NODE MODIFIERS **/

node DHBGraph::addNode() {
    m_dhb_graph.resize(m_dhb_graph.vertices_count() + 1);
    return m_dhb_graph.vertices_count() - 1;
}

node DHBGraph::addNodes(count numberOfNewNodes) {
    m_dhb_graph.resize(m_dhb_graph.vertices_count() + numberOfNewNodes);
    return m_dhb_graph.vertices_count() - 1;
}

void DHBGraph::removeNode(node v) {
    assert(v < m_dhb_graph.vertices_count());
    if (directed) {
        std::vector<std::pair<node, node>> edgesToRemove;
        forEdges([&](node from, node to) {
            if (from == v || to == v) {
                edgesToRemove.emplace_back(from, to);
            }
        });

        for (auto const &edge : edgesToRemove) {
            removeEdge(edge.first, edge.second);
        }
    } else {
        std::vector<node> neighbors;
        forNeighborsOf(v, [&](node u) { neighbors.push_back(u); });
        for (auto const &neighbor : neighbors) {
            removeEdge(v, neighbor);
            removeEdge(neighbor, v);
        }
    }
}

/** NODE PROPERTIES **/

edgeweight DHBGraph::weightedDegree(node u, bool countSelfLoopsTwice) const {
    return computeWeightedDegree(u, false, countSelfLoopsTwice);
}

edgeweight DHBGraph::weightedDegreeIn(node u, bool countSelfLoopsTwice) const {
    return computeWeightedDegree(u, true, countSelfLoopsTwice);
}

bool DHBGraph::hasNode(node v) const noexcept {
    return (v < m_dhb_graph.vertices_count());
}

/** EDGE MODIFIERS **/

// This interface is private. We are hiding the edge id as a parameter from the
// user since we're managing edge IDs internally.
bool DHBGraph::addEdge(node const u, node const v, edgeweight const ew, edgeid const id,
                       bool do_update) {
    assert(u < m_dhb_graph.vertices_count());
    assert(v < m_dhb_graph.vertices_count());

    auto const edge_weight = weighted ? ew : defaultEdgeWeight;

    // TODO: Why is dhb::EdgeID of unsigned int? That would be 32 bit unsigned
    // integer.. why is it not a 64 bit unsigned integer.
    //
    // For now, we will accept the "narrowing conversion of ‘id’ from ‘const
    // edgeid’ {aka ‘const long unsigned int’} to ‘dhb::EdgeID’ {aka ‘unsigned
    // int’}".
    EdgeData const edge_data{edge_weight, static_cast<dhb::EdgeID>(id)};

    bool is_inserted_in = false;
    bool is_inserted_out = false;
    bool is_inserted = false;

    // Insert edges depending on whether the graph is directed or undirected
    if (!directed) {
        is_inserted_out = m_dhb_graph.insert(u, v, edge_data, do_update);
        is_inserted_in = m_dhb_graph.insert(v, u, edge_data, do_update);
    } else {
        is_inserted = m_dhb_graph.insert(u, v, edge_data, do_update);
    }

    bool const undirected_insertion_success = is_inserted_in && is_inserted_out && !directed;
    bool const directed_insertion_success = is_inserted && directed;

    if (undirected_insertion_success || directed_insertion_success) {
        if (u == v) {
            ++storedNumberOfSelfLoops;
        }
        if (edgesIndexed) {
            omega++;
        }
        m++;
        return true;
    }
    return false;
}

// This is the public interface.
bool DHBGraph::addEdge(node const u, node const v, edgeweight const ew, bool do_update) {
    assert(u < m_dhb_graph.vertices_count());
    assert(v < m_dhb_graph.vertices_count());

    edgeid const id = edgesIndexed ? omega : 0;
    auto const applied_edge_weight = weighted ? ew : defaultEdgeWeight;
    return addEdge(u, v, applied_edge_weight, id, do_update);
}

void DHBGraph::addEdges(std::vector<WeightedEdge> &&weighted_edges, bool const do_update) {
    auto cmp = [](WeightedEdge const &a, WeightedEdge const &b) { return a.u < b.u; };
    auto get_source_f = [](WeightedEdge const &e) { return e.u; };

    auto insert_edge_f = [&](WeightedEdge const &e) {
        EdgeData data{e.weight, 0};
        auto [_, inserted] = m_dhb_graph.neighbors(e.u).insert(e.v, data);
        if (inserted) {
            storedNumberOfSelfLoops += uint64_t(e.u == e.v);
        }
    };

    auto insert_and_update_f = [&](WeightedEdge const &e) {
        EdgeData data{e.weight, 0};
        auto [it, inserted] = m_dhb_graph.neighbors(e.u).insert(e.v, data);
        if (inserted) {
            storedNumberOfSelfLoops += uint64_t(e.u == e.v);
        } else {
            it->data() = EdgeData{e.weight, it->data().id};
        }
    };

    dhb::BatchParallelizer<WeightedEdge> par;

    // TODO: we should not differ between updating and not updating edges within
    // the parallel for loop. This may slow down our inserts.
    auto processEdge = [&](WeightedEdge const &e) {
        if (do_update) {
            insert_and_update_f(e);
        } else {
            insert_edge_f(e);
        }
    };

    par(std::begin(weighted_edges), std::end(weighted_edges), get_source_f, cmp, processEdge);

    if (!isDirected()) {
        std::vector<WeightedEdge> weighted_edges_v_to_u(weighted_edges.size());

#pragma omp parallel for schedule(guided)
        for (omp_index i = 0u; i < static_cast<omp_index>(weighted_edges_v_to_u.size()); ++i) {
            weighted_edges_v_to_u[i].u = weighted_edges[i].v;
            weighted_edges_v_to_u[i].v = weighted_edges[i].u;
            weighted_edges_v_to_u[i].weight = weighted_edges[i].weight;
        }

        par(std::begin(weighted_edges_v_to_u), std::end(weighted_edges_v_to_u), get_source_f, cmp,
            processEdge);
    }
}

void DHBGraph::addEdges(std::vector<Edge> &&edges, bool do_update) {
    std::vector<WeightedEdge> weighted_edges;
    weighted_edges.reserve(edges.size());
    for (auto &edge : edges) {
        WeightedEdge w_edge = WeightedEdge{edge.u, edge.v, defaultEdgeWeight};
        weighted_edges.push_back(w_edge);
    }

    addEdges(std::move(weighted_edges), do_update);
}

template <typename T>
void erase(node u, index idx, std::vector<std::vector<T>> &vec) {
    vec[u][idx] = vec[u].back();
    vec[u].pop_back();
}

bool DHBGraph::removeEdge(node u, node v) {
    bool const is_loop = (u == v);
    bool is_removed_in = false;
    bool is_removed_out = false;
    bool is_removed = false;

    if (!directed) {
        is_removed_out = m_dhb_graph.removeEdge(u, v);
        is_removed_in = m_dhb_graph.removeEdge(v, u);
    } else {
        is_removed = m_dhb_graph.removeEdge(u, v);
    }

    bool const undirected_removal_success = (is_removed_in || is_removed_out) && !isDirected();
    bool const directed_removal_success = is_removed && isDirected();

    if (undirected_removal_success || directed_removal_success) {
        if (is_loop) {
            storedNumberOfSelfLoops--;
        }
        m--;
        return true;
    }
    return false;
}

void DHBGraph::removeAllEdges() {
    m_dhb_graph = dhb::Matrix<EdgeData>(m_dhb_graph.vertices_count());
    m = 0;
    storedNumberOfSelfLoops = 0;
    omega = 0;
}

void DHBGraph::removeSelfLoops() {
    parallelForNodes([&](node const u) {
        auto isSelfLoop = [u](node const v) { return u == v; };
        removeAdjacentEdges(u, isSelfLoop);
        if (isDirected()) {
            removeAdjacentEdges(u, isSelfLoop, true);
        }
    });

    m -= storedNumberOfSelfLoops;
    storedNumberOfSelfLoops = 0;
}

count DHBGraph::numberOfEdges() const noexcept {
    return m;
}

index DHBGraph::upperNodeIdBound() const noexcept {
    return m_dhb_graph.vertices_count();
}

count DHBGraph::numberOfNodes() const noexcept {
    return m_dhb_graph.vertices_count();
}

count DHBGraph::degree(node v) const {
    assert(v < m_dhb_graph.vertices_count());
    return m_dhb_graph.degree(v);
}

count DHBGraph::degreeIn(node v) const {
    assert(v < m_dhb_graph.vertices_count());
    if (directed) {
        count in_degree = 0;
        parallelForEdges([&](node, node to) {
            if (v == to) {
#pragma omp atomic
                ++in_degree;
            }
        });
        return in_degree;
    }
    return degreeOut(v);
}

void DHBGraph::swapEdge(node source_a, node target_a, node source_b, node target_b) {
    assert(source_a < m_dhb_graph.vertices_count());
    assert(target_a < m_dhb_graph.vertices_count());
    assert(source_b < m_dhb_graph.vertices_count());

    auto neigh_a = m_dhb_graph.neighbors(source_a);
    auto edge_a = neigh_a.iterator_to(target_a);
    if (edge_a == neigh_a.end()) {
        throw std::runtime_error("Edge (source_a, target_a) does not exist!");
    }

    auto neigh_b = m_dhb_graph.neighbors(source_b);
    auto edge_b = neigh_b.iterator_to(target_b);
    if (edge_b == neigh_b.end()) {
        throw std::runtime_error("Edge (source_b, target_b) does not exist!");
    }

    if (neigh_a.exists(target_b)) {
        throw std::runtime_error("Edge (source_a, target_b) does already exist!");
    }

    if (neigh_b.exists(target_a)) {
        throw std::runtime_error("Edge (source_b, target_a) does already exist!");
    }

    addEdge(source_a, target_b, edge_a->data().weight, edge_a->data().id);
    addEdge(source_b, target_a, edge_b->data().weight, edge_b->data().id);

    removeEdge(source_a, target_a);
    removeEdge(source_b, target_b);
}

bool DHBGraph::hasEdge(node u, node v) const noexcept {
    if (u < m_dhb_graph.vertices_count()) {
        return m_dhb_graph.neighbors(u).exists(v);
    }

    return false;
}

/** EDGE ATTRIBUTES **/

edgeweight DHBGraph::weight(node u, node v) const {
    auto it = m_dhb_graph.neighbors(u).iterator_to(v);
    if (it != m_dhb_graph.neighbors(u).end()) {
        return isWeighted() ? it->data().weight : defaultEdgeWeight;
    }
    return nullWeight;
}

void DHBGraph::setWeight(node u, node v, edgeweight ew) {
    if (!hasEdge(u, v)) {
        // edge does not exist, create it, but warn user
        TRACE("Setting edge weight of a nonexisting edge will create the edge.");
        addEdge(u, v, ew);
        return;
    }

    m_dhb_graph.neighbors(u).iterator_to(v)->data().weight = isWeighted() ? ew : defaultEdgeWeight;
    if (!directed) {
        m_dhb_graph.neighbors(v).iterator_to(u)->data().weight =
            isWeighted() ? ew : defaultEdgeWeight;
    }
}

void DHBGraph::increaseWeight(node u, node v, edgeweight ew) {
    if (!weighted) {
        throw std::runtime_error("Cannot increase edge weight in unweighted graph.");
    }

    if (!hasEdge(u, v)) {
        // edge does not exist, create it, but warn user
        TRACE("Increasing the weight of a nonexistent edge will create the edge with the proposed "
              "increased weight.");
        addEdge(u, v, ew);
        return;
    }
    edgeweight const new_ew = ew + weight(u, v);
    m_dhb_graph.neighbors(u).iterator_to(v)->data().weight = new_ew;
    if (!directed) {
        m_dhb_graph.neighbors(v).iterator_to(u)->data().weight = new_ew;
    }
}

void DHBGraph::setWeightAtIthNeighbor(Unsafe, node u, index i, edgeweight ew) {
    m_dhb_graph.neighbors(u)[i].data().weight = ew;
}

index DHBGraph::indexOfNeighbor(node u, node v) const {
    for (index i = 0; i < m_dhb_graph.neighbors(u).degree(); i++) {
        node current = m_dhb_graph.neighbors(u)[i].vertex();
        if (v == current) {
            return i;
        }
    }
    return none;
}

node DHBGraph::getIthNeighbor(node u, index i) const {
    if (!ithNeighborExists(u, i)) {
        return none;
    }
    return m_dhb_graph.neighbors(u)[i].vertex();
}

edgeweight DHBGraph::getIthNeighborWeight(node u, index i) const {
    if (!ithNeighborExists(u, i)) {
        return defaultEdgeWeight;
    }
    return isWeighted() ? m_dhb_graph.neighbors(u)[i].data().weight : defaultEdgeWeight;
}

std::pair<node, edgeweight> DHBGraph::getIthNeighborWithWeight(node u, index i) const {
    if (!ithNeighborExists(u, i)) {
        return {none, defaultEdgeWeight};
    }
    node const neighbor = m_dhb_graph.neighbors(u)[i].vertex();
    double const weight =
        isWeighted() ? m_dhb_graph.neighbors(u)[i].data().weight : defaultEdgeWeight;
    return {neighbor, weight};
}

std::pair<node, edgeid> DHBGraph::getIthNeighborWithId(node u, index i) const {
    assert(hasEdgeIds());

    if (!ithNeighborExists(u, i)) {
        return {none, none};
    }
    node const neighbor = m_dhb_graph.neighbors(u)[i].vertex();
    edgeid const id = m_dhb_graph.neighbors(u)[i].data().id;
    return {neighbor, id};
}

/** SUMS **/

edgeweight DHBGraph::totalEdgeWeight() const noexcept {
    if (weighted)
        return sumForEdgesParallel([](node, node, edgeweight ew) { return ew; });
    else
        return numberOfEdges() * defaultEdgeWeight;
}

} // namespace NetworKit

#include <networkit/matching/BMatching.hpp>

namespace NetworKit {

BMatching::BMatching(const Graph &G, const std::vector<count> &b)
    : G(&G), b(b), matches(G.numberOfNodes()) {}

bool BMatching::isProper() const {
    // check if entries are symmetric and every pair exists as an edge
    for (node v : G->nodeRange()) {
        for (node w : matches[v]) {
            if (matches[w].find(v) == matches[w].end()) {
                DEBUG("node ", v, " is not symmetrically matched");
                return false;
            }
            if ((v != w) && !G->hasEdge(v, w)) {
                DEBUG("matched pair (", v, ",", w, ") is not an edge");
                return false;
            }
        }
    }
    return true;
}

void BMatching::match(node u, node v) {
    assert(matches[u].size() <= b[u]);
    assert(matches[v].size() <= b[v]);
    matches[u].insert(v);
    matches[v].insert(u);
}

void BMatching::unmatch(node u, node v) {
    matches[u].erase(v);
    matches[v].erase(u);
}

bool BMatching::isUnmatched(node u) const {
    return matches[u].empty();
}

bool BMatching::areMatched(node u, node v) const {
    return matches[u].find(v) != matches[u].end();
}

count BMatching::size() const {
    double size = G->parallelSumForNodes([&](node v) { return !isUnmatched(v); });
    return static_cast<count>(size / 2);
}

edgeweight BMatching::weight() const {
    return G->parallelSumForNodes([&](node v) {
        edgeweight weight_per_node = 0.0;
        for (node u : matches[v]) {
            if (v < u) {
                weight_per_node += G->weight(v, u);
            }
        }
        return weight_per_node;
    });
}

const std::vector<std::unordered_set<node>> &BMatching::getMatches() const {
    return matches;
}

const std::vector<count> &BMatching::getB() const {
    return b;
}

void BMatching::reset() {
    for (auto &nodeMatches : matches) {
        nodeMatches.clear();
    }
}

} // namespace NetworKit

#include <networkit/matching/BMatching.hpp>

namespace NetworKit {

BMatching::BMatching(const std::vector<count> &b, count z) : b(b) {
    matches.resize(z);
}

bool BMatching::isProper(const Graph &G) const {
    // check if entries are symmetric and every pair exists as an edge
    for (node v : G.nodeRange()) {
        for (node w : matches.at(v)) {
            if (matches.at(w).find(v) == matches.at(w).end()) {
                DEBUG("node ", v, " is not symmetrically matched");
                return false;
            }
            if ((v != w) && !G.hasEdge(v, w)) {
                DEBUG("matched pair (", v, ",", w, ") is not an edge");
                return false;
            }
        }
    }
    return true;
}

void BMatching::match(node u, node v) {
    assert(matches.at(u).size() <= b.at(u));
    assert(matches.at(v).size() <= b.at(v));
    matches.at(u).insert(v);
    matches.at(v).insert(u);
}

void BMatching::unmatch(node u, node v) {
    matches.at(u).erase(v);
    matches.at(v).erase(u);
}

bool BMatching::isUnmatched(node u) const {
    return matches.at(u).empty();
}

bool BMatching::areMatched(node u, node v) const {
    return matches.at(u).find(v) != matches.at(u).end();
}

count BMatching::size(const Graph &G) const {
    count size = 0;
    G.forNodes([&](node v) {
        if (!isUnmatched(v)) {
            ++size;
        }
    });
    return size / 2;
}

edgeweight BMatching::weight(const Graph &G) const {
    return G.parallelSumForNodes([&](node v) {
        edgeweight weight_per_node = 0.0;
        for (auto u : matches.at(v)) {
            if (v < u) {
                weight_per_node += G.weight(v, u);
            }
        }
        return weight_per_node;
    });
}

const std::vector<std::set<node>> &BMatching::getMatches() const {
    return matches;
}

std::vector<count> BMatching::getB() const {
    return b;
}

} // namespace NetworKit

#include <networkit/matching/BSuitorMatcher.hpp>

#include <algorithm>
#include <iomanip>
#include <ranges>
#include <stdexcept>

namespace NetworKit {
BSuitorMatcher::BSuitorMatcher(const Graph &G, const std::vector<count> &b) : BMatcher(G, b), b(b) {
    if (G.numberOfSelfLoops() > 0)
        throw std::runtime_error("This algorithm does not support graphs with self-loops.");

    if (G.isDirected())
        throw std::runtime_error("This algorithm does not support directed graphs.");

    if (b.size() != G.numberOfNodes())
        throw std::runtime_error(
            "The number of b values does not match the number of nodes in this graph.");

    if (!G.isWeighted())
        throw std::runtime_error("This algorithm does not support unweighted graphs, use "
                                 "GraphTools::toWeighted to convert the graph.");

    const auto n = G.upperNodeIdBound();

    if (n != G.numberOfNodes())
        throw std::runtime_error(
            "The graph needs to be compact in order to calculate the b-matching, e.g., preprocess "
            "the graph with GraphTools::getCompactedGraph.");

    suitors.reserve(n);
    proposed.reserve(n);
    for (index i = 0; i < n; ++i) {
        suitors.emplace_back(b[i]);
        proposed.emplace_back(b[i]);
    }
}

BSuitorMatcher::BSuitorMatcher(const Graph &G, count b)
    : BSuitorMatcher(G, std::vector<count>(G.numberOfNodes(), b)) {}

void BSuitorMatcher::findSuitors(node cur) {
    for (index i = 0; i < b[cur]; i++) {
        auto [pref, heaviest] = findPreferred(cur);
        if (pref != none) {
            makeSuitor(cur, heaviest, pref);
        }
    }
}

BSuitorMatcher::MatchingNode BSuitorMatcher::findPreferred(node u) {
    MatchingNode best = MatchingNode{none, 0};

    auto hasProposedTo = [&](node x) -> bool {
        return std::ranges::any_of(proposed[u].partners,
                                   [x](const MatchingNode &y) { return y.id == x; });
    };

    for (auto [v, weight] : G->weightNeighborRange(u)) {
        const MatchingNode w = MatchingNode(v, weight);
        if (hasProposedTo(w.id))
            continue;
        if (w > best) {
            const edgeweight n_suitor_weight = suitors[w.id].min.weight;

            if (w.weight > n_suitor_weight
                || (w.weight == n_suitor_weight && u < suitors[w.id].min.id)) {
                best = w;
            }
        }
    }
    return best;
}

void BSuitorMatcher::makeSuitor(node u, edgeweight w, node v) {
    auto smallest = suitors[v].insert(MatchingNode(u, w));
    proposed[u].insert(MatchingNode(v, w));

    if (smallest.id != none) {
        proposed[smallest.id].remove(v);
        auto [pref, heaviest] = findPreferred(smallest.id);
        if (pref != none) {
            makeSuitor(smallest.id, heaviest, pref);
        }
    }
}

bool BSuitorMatcher::isSymmetrical() const {
    bool sym = true;
    auto matchedSymmetrical = [&](node x, node y) -> bool {
        return suitors[x].hasPartner(y) == suitors[y].hasPartner(x);
    };

    for (node u : G->nodeRange()) {
        for (node v = u; v < G->upperNodeIdBound(); ++v) {
            if (!matchedSymmetrical(u, v)) {
                sym = false;
                break;
            }
        }
    }

    return sym;
}

void BSuitorMatcher::buildBMatching() {
    bMatch.reset();
    G->forNodes([&](node x) {
        assert(suitors[x].partners.size() <= b.at(x));
        for (MatchingNode y : suitors[x].partners) {
            if (y.id != none && x < y.id) {
                bMatch.match(x, y.id);
            }
        }
    });
}

void BSuitorMatcher::run() {
    G->forNodes([&](node u) { findSuitors(u); });
    buildBMatching();
    hasRun = true;
}
} // namespace NetworKit

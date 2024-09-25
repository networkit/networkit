#include <networkit/matching/BSuitorMatcher.hpp>

#include <algorithm>
#include <iomanip>
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

    const auto n = G.upperNodeIdBound();

    if (n != G.numberOfNodes())
        throw std::runtime_error(
            "The graph needs to be compact in order to calculate the b-matching, e.g., preprocess the graph with GraphTools::getCompactedGraph.");

    suitors.reserve(n);
    proposed.reserve(n);
    for (index i = 0; i < n; ++i) {
        suitors.emplace_back(std::make_unique<MatchingNodeInfo>(b.at(i)));
        proposed.emplace_back(std::make_unique<MatchingNodeInfo>(b.at(i)));
    }
}

BSuitorMatcher::BSuitorMatcher(const Graph &G, count b)
    : BSuitorMatcher(G, std::vector<count>(G.numberOfNodes(), b)) {}

BSuitorMatcher::BSuitorMatcher(const Graph &G, std::string_view &path)
    : BSuitorMatcher(G, readBValuesFromFile(G.numberOfNodes(), path)) {}

std::vector<count> BSuitorMatcher::readBValuesFromFile(count size, std::string_view &path) const {
    std::vector<count> b;
    b.reserve(size);
    std::ifstream file(path.data());
    std::string line;
    int line_number = 1;

    while (std::getline(file, line)) {
        std::istringstream istring(line);
        int val;
        if (!(istring >> val)) {
            throw std::runtime_error("File " + std::string{path}
                                     + " contains an invalid value in line "
                                     + std::to_string(line_number) + ".");
        }
        if (istring >> val) {
            throw std::runtime_error("File " + std::string{path}
                                     + " contains multiple values in line "
                                     + std::to_string(line_number) + ".");
        }
        if (val < 0) {
            throw std::runtime_error("File " + std::string{path}
                                     + " contains a negative value in line "
                                     + std::to_string(line_number) + ".");
        }
        b.emplace_back(val);
        ++line_number;
    }
    if (b.size() != size) {
        throw std::runtime_error("The number of values in file " + std::string{path}
                                 + " does not match the number of nodes in this graph.");
    }
    return b;
}

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
        return std::any_of(proposed[u]->partners.begin(), proposed[u]->partners.end(),
                           [x](const MatchingNode &y) { return y.id == x; });
    };

    for (auto [v, weight] : G->weightNeighborRange(u)) {
        const MatchingNode w = MatchingNode(v, weight);
        if (hasProposedTo(w.id))
            continue;
        if (w > best) {
            const edgeweight n_suitor_weight = suitors[w.id]->min.weight;

            if (w.weight > n_suitor_weight
                || (w.weight == n_suitor_weight && u < suitors[w.id]->min.id)) {
                best = w;
            }
        }
    }
    return best;
}

void BSuitorMatcher::makeSuitor(node u, edgeweight w, node v) {
    auto smallest = suitors[v]->popMinIfFull();
    suitors[v]->insert(MatchingNode(u, w));
    proposed[u]->insert(MatchingNode(v, w));

    if (smallest.id != none) {
        proposed[smallest.id]->remove(v);
        auto [pref, heaviest] = findPreferred(smallest.id);
        if (pref != none) {
            makeSuitor(smallest.id, heaviest, pref);
        }
    }
}

bool BSuitorMatcher::isSymmetrical() const {
    bool sym = true;
    auto matchedSymmetrical = [&](node x, node y) -> bool {
        return suitors[x]->hasPartner(y) == suitors[y]->hasPartner(x);
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
    G->forNodes([&](node x) {
        assert(suitors[x]->partners.size() <= b.at(x));
        for (MatchingNode y : suitors[x]->partners) {
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

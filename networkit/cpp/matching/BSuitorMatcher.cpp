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
    suitors.reserve(n);
    proposed.reserve(n);
    for (index i = 0; i < n; i++) {
        suitors.emplace_back(std::make_unique<MatchingNodeInfo>(b.at(i)));
        proposed.emplace_back(std::make_unique<MatchingNodeInfo>(b.at(i)));
    }
}

BSuitorMatcher::BSuitorMatcher(const Graph &G, count b)
    : BSuitorMatcher(G, std::vector<count>(G.numberOfNodes(), b)) {}

BSuitorMatcher::BSuitorMatcher(const Graph &G, const std::string &path)
    : BSuitorMatcher(G, readBValuesFromFile(G.numberOfNodes(), path)) {}

std::vector<count> BSuitorMatcher::readBValuesFromFile(count size, const std::string &path) const {
    std::vector<count> b;
    b.reserve(size);
    std::ifstream file(path);
    std::string line;
    int line_number = 1;

    while (std::getline(file, line)) {
        std::istringstream istring(line);
        int val;
        if (!(istring >> val)) {
            throw std::runtime_error("File " + path + " contains an invalid value in line "
                                     + std::to_string(line_number) + ".");
        }
        if (istring >> val) {
            throw std::runtime_error("File " + path + " contains multiple values in line "
                                     + std::to_string(line_number) + ".");
        }
        if (val < 0) {
            throw std::runtime_error("File " + path + " contains a negative value in line "
                                     + std::to_string(line_number) + ".");
        }
        b.emplace_back(val);
        line_number++;
    }
    if (b.size() != size) {
        throw std::runtime_error("The number of values in file " + path
                                 + " does not match the number of nodes in this graph.");
    }
    return b;
}

void BSuitorMatcher::findSuitors(node cur) {
    for (index i = 0; i < b.at(cur); i++) {
        auto [pref, heaviest] = findPreferred(cur);
        if (pref != none) {
            makeSuitor(cur, heaviest, pref);
        }
    }
}

BSuitorMatcher::MatchingNode BSuitorMatcher::findPreferred(node u) {
    MatchingNode best = MatchingNode{none, 0};

    auto hasProposedTo = [&](node x) -> bool {
        return std::any_of(proposed.at(u)->partners.begin(), proposed.at(u)->partners.end(),
                           [x](const MatchingNode &y) { return y.id == x; });
    };

    for (auto n : G->weightNeighborRange(u)) {
        const MatchingNode v = MatchingNode(n.first, n.second);
        if (!hasProposedTo(v.id)) {
            if (v.weight > best.weight || (v.weight == best.weight && v.id < best.id)) {
                const auto n_suitor_weight = suitors.at(v.id)->min.weight;

                if (v.weight > n_suitor_weight
                    || (v.weight == n_suitor_weight && u < suitors.at(v.id)->min.id)) {
                    best = v;
                }
            }
        }
    }
    return best;
}

void BSuitorMatcher::makeSuitor(node u, edgeweight w, node v) {
    auto smallest = suitors.at(v)->popMinIfFull();
    suitors.at(v)->insert(MatchingNode(u, w));
    proposed.at(u)->insert(MatchingNode(v, w));

    if (smallest.id != none) {
        proposed.at(smallest.id)->remove(v);
        auto [pref, heaviest] = findPreferred(smallest.id);
        if (pref != none) {
            makeSuitor(smallest.id, heaviest, pref);
        }
    }
}

bool BSuitorMatcher::isSymmetrical() const {
    bool sym = true;
    auto matchedSymmetrical = [&](node x, node y) -> bool {
        return suitors.at(x)->hasPartner(y) == suitors.at(y)->hasPartner(x);
    };

    G->forNodes([&](node u) {
        G->forNodes([&](node v) {
            if (u > v && !matchedSymmetrical(u, v)) {
                sym = false;
            }
        });
    });
    return sym;
}

void BSuitorMatcher::buildBMatching() {
    if (!hasRun) {
        throw std::runtime_error("Call run() before creating b-matching.");
    }
    G->forNodes([&](node x) {
        assert(suitors.at(x)->partners.size() <= b.at(x));
        for (auto y : suitors.at(x)->partners) {
            if (y.id != none && x < y.id) {
                M.match(x, y.id);
            }
        }
    });
}

void BSuitorMatcher::run() {
    G->forNodes([&](node u) { findSuitors(u); });
    hasRun = true;
}
} // namespace NetworKit

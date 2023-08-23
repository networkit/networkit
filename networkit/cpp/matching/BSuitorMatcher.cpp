#include <networkit/matching/BSuitorMatcher.hpp>

#include <algorithm>
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
}

BSuitorMatcher::BSuitorMatcher(const Graph &G, count b)
    : BSuitorMatcher(G, std::vector<count>(G.numberOfNodes(), b)) {}

BSuitorMatcher::BSuitorMatcher(const Graph &G, const std::string &path)
    : BSuitorMatcher(G, readBValuesFromFile(G.numberOfNodes(), path)) {}

void BSuitorMatcher::findSuitors(node u) {
    for (count i = 0; i < b.at(u); i++) {
        auto x = findPreffered(u);
        if (x != none)
            makeSuitor(u, x);
    }
}

node BSuitorMatcher::findPreffered(node y) {
    node x = none;
    edgeweight heaviest = 0;

    auto hasProposedTo = [&](node n) -> bool {
        return (std::find(proposed[y].begin(), proposed[y].end(), n) != proposed[y].end());
    };

    for (auto v : G->weightNeighborRange(y)) {
        if (!hasProposedTo(v.first)) {
            const edgeweight weight = v.second;
            if (weight > heaviest || (weight == heaviest && v.first < x)) {
                if (weight > G->weight(v.first, suitors[v.first].back())
                    || (weight == G->weight(v.first, suitors[v.first].back())
                        && y < suitors[v.first].back())) {
                    x = v.first;
                    heaviest = weight;
                }
            }
        }
    }
    return x;
}

void BSuitorMatcher::makeSuitor(node u, node x) {
    auto y = suitors[x].back();
    sortInsert(suitors[x], x, u);
    auto i = findFirstFreeIndex(proposed[u]);
    proposed[u][i] = x;

    if (y != none) {
        sortRemove(proposed[y], x);
        auto z = findPreffered(y);

        if (z != none) {
            makeSuitor(y, z);
        }
    }
}

std::vector<count> BSuitorMatcher::readBValuesFromFile(count size, const std::string &path) const {
    std::vector<count> b;
    b.reserve(size);

    std::ifstream file(path);
    std::string line;
    int val;

    while (std::getline(file, line)) {
        std::istringstream istring(line);

        if (!(istring >> val) || val < 0) {
            throw std::runtime_error("Invalid value (" + std::to_string(val) + ") found in file "
                                     + path + ".");
        }
        b.push_back(val);
    }

    if (b.size() != size) {
        throw std::runtime_error("The number of values in file " + path
                                 + " does not match the number of nodes in this graph.");
    }

    return b;
}

void BSuitorMatcher::sortInsert(std::vector<node> &nodes, node u, node v) {
    auto i = std::find_if(nodes.begin(), nodes.end(), [&](const node &y) {
        return (G->weight(u, v) > G->weight(u, y))
               || ((G->weight(u, v) == G->weight(u, y)) && (v < y));
    });
    nodes.insert(i, v);
    nodes.pop_back();
}

void BSuitorMatcher::sortRemove(std::vector<node> &nodes, const node u) {
    nodes.erase(std::remove(nodes.begin(), nodes.end(), u), nodes.end());
    nodes.push_back(none);
}

index BSuitorMatcher::findFirstFreeIndex(const std::vector<node> &nodes) const {
    return findIndexOf(nodes, none);
}

index BSuitorMatcher::findIndexOf(const std::vector<node> &nodes, node x) const {
    auto i = std::find(nodes.begin(), nodes.end(), x);
    return i == nodes.end() ? none : i - nodes.begin();
}

void BSuitorMatcher::checkSymmetry() const {
    auto areMatchedSymmetrical = [&](node u, node v) -> void {
        [[maybe_unused]] auto i_1 = std::find(this->suitors[u].begin(), this->suitors[u].end(), v)
                                    != this->suitors[u].end();
        [[maybe_unused]] auto i_2 = std::find(this->suitors[v].begin(), this->suitors[v].end(), u)
                                    != this->suitors[v].end();
        assert(i_1 == i_2);
    };

    G->forNodes([&](node u) { G->forNodes([&](node v) { areMatchedSymmetrical(u, v); }); });
}

void BSuitorMatcher::run() {
    const auto n = G->upperNodeIdBound();
    for (index i = 0; i < n; i++) {
        auto v = std::vector<node>(b.at(i), none);
        suitors.emplace_back(v);
        proposed.emplace_back(v);
    }

    G->forNodes([&](node u) { findSuitors(u); });
    checkSymmetry();

    // TODO make parallel
    G->forNodes([&suitors = suitors, &M = M](node u) {
        for (auto p : suitors[u]) {
            if (p != none && u < p) { // Ensure we match a pair of nodes only once
                M.match(u, p);
            }
        }
    });
    hasRun = true;
}
} // namespace NetworKit
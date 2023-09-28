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
}

BSuitorMatcher::BSuitorMatcher(const Graph &G, count b)
    : BSuitorMatcher(G, std::vector<count>(G.numberOfNodes(), b)) {}

BSuitorMatcher::BSuitorMatcher(const Graph &G, const std::string &path)
    : BSuitorMatcher(G, readBValuesFromFile(G.numberOfNodes(), path)) {}

void BSuitorMatcher::findSuitors(node u) {
    for (index i = 0; i < b.at(u); i++) {
        auto [x, w] = findPreferred(u);
        if (x != none) {
            makeSuitor(u, w, x);
        }
    }
}

Node BSuitorMatcher::findPreferred(node y) {
    Node best = Node{none, 0};

    auto hasProposedTo = [&](node n) -> bool {
        return std::any_of(T.at(y)->list.begin(), T.at(y)->list.end(),
                           [n](const Node &p) { return p.id == n; });
    };

    for (auto v : G->weightNeighborRange(y)) {
        const Node n = Node(v.first, v.second);
        if (!hasProposedTo(n.id)) {
            if (n.weight > best.weight || (n.weight == best.weight && n.id < best.id)) {
                const auto n_suitor_weight = S.at(n.id)->min.weight;

                if (n.weight > n_suitor_weight
                    || (n.weight == n_suitor_weight && y < S.at(n.id)->min.id)) {
                    best = n;
                }
            }
        }
    }
    return best;
}

void BSuitorMatcher::makeSuitor(node u, edgeweight w, node x) {
    auto y = S.at(x)->popMinIfFull();
    S.at(x)->insert(Node(u, w));
    T.at(u)->insert(Node(x, w));

    if (y != none) {
        T.at(y)->remove(Node(x, w));
        auto [z, w] = findPreferred(y);
        if (z != none) {
            makeSuitor(y, w, z);
        }
    }
}

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
        b.push_back(val);
        line_number++;
    }
    if (b.size() != size) {
        throw std::runtime_error("The number of values in file " + path
                                 + " does not match the number of nodes in this graph.");
    }
    return b;
}

bool BSuitorMatcher::isSymmetrical() const {
    bool sym = true;
    auto areMatchedSymmetrical = [&](node u, node v) -> bool {
        auto i_1 = std::find_if(S.at(u)->list.begin(), S.at(u)->list.end(),
                                [v](const Node &p) { return p.id == v; })
                   != S.at(u)->list.end();
        auto i_2 = std::find_if(S.at(v)->list.begin(), S.at(v)->list.end(),
                                [u](const Node &p) { return p.id == u; })
                   != S.at(v)->list.end();
        return (i_1 == i_2);
    };

    G->forNodes([&](node u) {
        G->forNodes([&](node v) {
            if (!areMatchedSymmetrical(u, v)) {
                sym = false;
            }
        });
    });
    return sym;
}

edgeweight BSuitorMatcher::getWeight() const {
    return w;
}

void BSuitorMatcher::run() {
    const auto n = G->upperNodeIdBound();
    S.reserve(n);
    T.reserve(n);
    for (index i = 0; i < n; i++) {
        Suitors *s = new Suitors(i, b.at(i));
        Proposed *p = new Proposed(i, b.at(i));
        S.emplace_back(s);
        T.emplace_back(p);
    }

    G->forNodes([&](node u) { findSuitors(u); });

#pragma omp parallel for reduction(+ : w)
    for (omp_index i = 0; i < static_cast<omp_index>(S.size()); i++) {
        w += S.at(i)->weight;
    }

    // TODO make parallel
    G->forNodes([&S = S, &M = M](node u) {
        for (auto v : S.at(u)->list) {
            if (v.id != none && u < v.id) { // Ensure we match a pair of nodes only once
                M.match(u, v.id);
            }
        }
    });

    for (auto s : S)
        delete s;
    for (auto t : T)
        delete t;

    hasRun = true;
}
} // namespace NetworKit
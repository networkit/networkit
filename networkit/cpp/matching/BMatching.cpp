#include <networkit/matching/BMatching.hpp>

namespace NetworKit {

BMatching::BMatching(const std::vector<count> &b, count z)
    : data(z, std::vector<node>(b.size(), none)), b(b) {}

bool BMatching::isProper(const Graph &G) const {
    // check if entries are symmetric
    for (node v : G.nodeRange()) {
        for (index i = 0; i < b.at(i); i++) {
            node w = data[v][i];
            if (w != none && std::find(data[w].begin(), data[w].end(), v) == data[w].end()) {
                DEBUG("node ", v, " is not symmetrically matched");
                return false;
            }
        }
    }

    // check if every pair exists as an edge
    for (node v : G.nodeRange()) {
        for (index i = 0; i < b.at(i); i++) {
            node w = data[v][i];
            if ((v != w) && (w != none) && !G.hasEdge(v, w)) {
                DEBUG("matched pair (", v, ",", w, ") is not an edge");
                return false;
            }
        }
    }

    return true;
}

void BMatching::match(node u, node v) {
    {
        auto i = findFirstFreeIndex(u);
        assert(i != none);
        data[u][i] = v;
    }
    {
        auto i = findFirstFreeIndex(v);
        assert(i != none);
        data[v][i] = u;
    }
}

void BMatching::unmatch(node u, node v) {
    {
        auto i = findIndexOf(u, v);
        assert(i != none);
        data[u][i] = none;
    }
    {
        auto i = findIndexOf(v, u);
        assert(i != none);
        data[v][i] = none;
    }
}

bool BMatching::isUnmatched(node u) const {
    return (this->data[u].front() == none);
}

bool BMatching::areMatched(node u, node v) const {
    return std::find(data[u].begin(), data[u].end(), v) != data[u].end();
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

std::vector<node> BMatching::mates(node v) const {
    return data[v];
}

edgeweight BMatching::weight(const Graph &G) const {
    edgeweight weight = 0;

    G.forNodes([&](node v) {
        if (!isUnmatched(v)) {
            for (auto u : mates(v)) {
                if (v < u)
                    weight += G.weight(v, u);
            }
        }
    });

    return weight;
}

std::vector<std::vector<node>> BMatching::getMatrix() const {
    return this->data;
}

std::vector<count> BMatching::getB() const {
    return b;
}

index BMatching::findFirstFreeIndex(node u) const {
    return findIndexOf(u, none);
}

index BMatching::findIndexOf(node u, node x) const {
    auto i = std::find(data[u].begin(), data[u].end(), x);
    return (i == std::end(data[u])) ? none : i - data[u].begin();
}

} // namespace NetworKit
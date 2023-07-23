#include <networkit/matching/BMatcher.hpp>

namespace NetworKit {

BMatcher::BMatcher(const Graph &G, int b) : G(&G), M(G.upperNodeIdBound(), b) {}

BMatching BMatcher::getBMatching() const {
    assureFinished();
    return M;
}

} /* namespace NetworKit */

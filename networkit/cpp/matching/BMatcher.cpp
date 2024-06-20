#include <networkit/matching/BMatcher.hpp>

namespace NetworKit {

BMatcher::BMatcher(const Graph &G, const std::vector<count> &b)
    : G(&G), M(b, G.upperNodeIdBound()) {}

BMatching BMatcher::getBMatching() const {
    assureFinished();
    return M;
}

} /* namespace NetworKit */

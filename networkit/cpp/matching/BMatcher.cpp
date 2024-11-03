#include <networkit/matching/BMatcher.hpp>

namespace NetworKit {

BMatcher::BMatcher(const Graph &G, const std::vector<count> &b) : G(&G), bMatch(G, b) {}

BMatching BMatcher::getBMatching() const {
    assureFinished();
    return bMatch;
}

} /* namespace NetworKit */

#ifndef NETWORKIT_MATCHING_B_MATCHER_HPP_
#define NETWORKIT_MATCHING_B_MATCHER_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/matching/BMatching.hpp>

namespace NetworKit {

/**
 * @ingroup matching
 * Base class for b-matching algorithms.
 */
class BMatcher : public Algorithm {
protected:
    const Graph *G;
    BMatching M;

public:
    /**
     * Constructs a new BMatcher.
     *
     * @param G
     * @param b
     */
    BMatcher(const Graph &G, const std::vector<count> &b);

    ~BMatcher() override = default;

    /**
     * Runs the b-matching algorithm on the stored graph.
     *
     */
    void run() override = 0;

    BMatching getBMatching() const;
};

} /* namespace NetworKit */
#endif // NETWORKIT_MATCHING_B_MATCHER_HPP_

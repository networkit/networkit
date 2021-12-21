/*
 * GraphletsCounter.hpp
 *
 * Created on: Dec 3, 2021
 *     Author: Bermudes
 */

#ifndef NETWORKIT_SUBGRAPHS_GRAPHLETS_COUNTER_HPP_
#define NETWORKIT_SUBGRAPHS_GRAPHLETS_COUNTER_HPP_

#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class GraphletsCounter final : public Algorithm {
public:
    GraphletsCounter(const Graph &G, unsigned K) noexcept;

    ~GraphletsCounter() noexcept override = default;

    void run() override;

    inline const std::vector<count> &getGraphletsCounts() const noexcept;

private:
    const Graph *G; // (style guide: keep pointers, not references)
    unsigned k;
    std::vector<count> counts;
};

const std::vector<count> &GraphletsCounter::getGraphletsCounts() const noexcept {
    return counts;
}

} // namespace NetworKit

#endif // NETWORKIT_SUBGRAPHS_GRAPHLETS_COUNTER_HPP_

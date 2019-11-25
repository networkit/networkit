/*
 * CurveballGlobalTradeGenerator.h
 *
 *  Created on: Jul 11, 2017
 *      Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */
// networkit-format

#ifndef NETWORKIT_RANDOMIZATION_CURVEBALL_GLOBAL_TRADE_GENERATOR_HPP_
#define NETWORKIT_RANDOMIZATION_CURVEBALL_GLOBAL_TRADE_GENERATOR_HPP_

#include <stdexcept>
#include <utility>
#include <vector>

#include <networkit/Globals.hpp>

namespace NetworKit {

class CurveballGlobalTradeGenerator {
public:
    using value_type = std::vector<std::pair<node, node>>;

protected:
    const node numNodes;
    const count numGlobalTrades;

public:
    CurveballGlobalTradeGenerator(const count run_length, const node num_nodes)
        : numNodes(num_nodes), numGlobalTrades(run_length) {
        if (num_nodes <= 1)
            throw std::runtime_error("At least two nodes are required.");
    }

    value_type generate() const;
};

} // namespace NetworKit

#endif // NETWORKIT_RANDOMIZATION_CURVEBALL_GLOBAL_TRADE_GENERATOR_HPP_

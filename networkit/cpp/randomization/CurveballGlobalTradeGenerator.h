/*
 * CurveballGlobalTradeGenerator.h
 *
 *  Created on: Jul 11, 2017
 *      Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */
#pragma once

#ifndef RANDOMIZATION_CURVEBALL_GLOBAL_TRADE_GENERATOR_H
#define RANDOMIZATION_CURVEBALL_GLOBAL_TRADE_GENERATOR_H

#include "../Globals.h"

namespace NetworKit {

class CurveballGlobalTradeGenerator {
public:
	using value_type = std::vector< std::pair<node, node> >;

protected:
	const node numNodes;
	const count numGlobalTrades;

public:
	CurveballGlobalTradeGenerator(const count run_length, const node num_nodes)
		: numNodes(num_nodes)
		, numGlobalTrades(run_length)
	{
		if (num_nodes <= 1)
			throw std::runtime_error("At least two nodes are required.");
	}

	value_type generate() const;
};

}

#endif // RANDOMIZATION_CURVEBALL_GLOBAL_TRADE_GENERATOR_H

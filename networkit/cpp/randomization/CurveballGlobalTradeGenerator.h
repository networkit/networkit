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
	const node _num_nodes;
	const count _run_length;
	const count _trades_per_run;

public:
	CurveballGlobalTradeGenerator(const count run_length, const node num_nodes)
		: _num_nodes(num_nodes)
		, _run_length(run_length)
		, _trades_per_run(static_cast<count>(num_nodes / 2))
	{
		if (num_nodes <= 1)
			throw std::runtime_error("At least two nodes are required.");
	}

	value_type generate() const;
};

}

#endif // RANDOMIZATION_CURVEBALL_GLOBAL_TRADE_GENERATOR_H

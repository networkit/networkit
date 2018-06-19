/*
 * CurveballUniformTradeGenerator.h
 *
 *  Created on: Jul 11, 2017
 *      Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */
#ifndef RANDOMIZATION_CURVEBALL_UNIFORM_TRADE_GENERATOR_H
#define RANDOMIZATION_CURVEBALL_UNIFORM_TRADE_GENERATOR_H

#include "../Globals.h"

namespace NetworKit {

class CurveballUniformTradeGenerator {
public:
	using value_type = std::vector<std::pair<node, node> >;

protected:
	const count _trade_num;
	const node _num_nodes;

public:
	CurveballUniformTradeGenerator(const count trade_num, const node num_nodes)
		: _trade_num(trade_num), _num_nodes(num_nodes)
	{
		if (num_nodes <= 1)
			throw std::runtime_error("At least two nodes are required.");
	}

	value_type generate() const;
};

}

#endif // RANDOMIZATION_CURVEBALL_UNIFORM_TRADE_GENERATOR_H

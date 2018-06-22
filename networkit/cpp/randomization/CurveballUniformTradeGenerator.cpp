/*
 * CurveballUniformTradeGenerator.cpp
 *
 *  Created on: Jul 11, 2017
 *  	Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */

#include <cassert>
#include <algorithm>
#include <vector>

#include "CurveballUniformTradeGenerator.h"

#include "../auxiliary/Random.h"


namespace NetworKit {

CurveballUniformTradeGenerator::value_type CurveballUniformTradeGenerator::generate() const {
	value_type trades_out;
	trades_out.clear();
	trades_out.reserve(numTrades);

	for (count t_id = 0; t_id < numTrades; t_id++) {
		const node fst = Aux::Random::integer(numNodes - 1);
		while (true) {
			const node snd = Aux::Random::integer(numNodes - 1);
			if (fst != snd) {
				trades_out.emplace_back(fst, snd);
				break;
			}
		}
	}

	return trades_out;
}

} // namespace NetworKit

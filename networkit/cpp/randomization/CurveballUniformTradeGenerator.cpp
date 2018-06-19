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
	value_type _trades_out;
	_trades_out.clear();
	_trades_out.reserve(_trade_num);

	for (count t_id = 0; t_id < _trade_num; t_id++) {
		const node fst = Aux::Random::integer(_num_nodes - 1);
		while (true) {
			const node snd = Aux::Random::integer(_num_nodes - 1);
			if (fst != snd) {
				_trades_out.emplace_back(fst, snd);
				break;
			}
		}
	}

	return _trades_out;
}

} // namespace NetworKit

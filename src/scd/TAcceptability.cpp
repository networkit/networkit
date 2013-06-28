/*
 * TAcceptability.cpp
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#include "TAcceptability.h"

namespace NetworKit {

TAcceptability::TAcceptability(const Graph& G, std::unordered_set<node>& community, std::unordered_set<node>& shell) : G(&G), community(&community), shell(&shell) {
	// TODO Auto-generated constructor stub

}

TAcceptability::~TAcceptability() {
	// TODO Auto-generated destructor stub
}

double TAcceptability::getValue(node v) {
	throw std::runtime_error("calling abstract method - needs to be redefined in subclass");
}

} /* namespace NetworKit */

NetworKit::TDummyAcceptability::TDummyAcceptability(const Graph& G,
		std::unordered_set<node>& community, std::unordered_set<node>& shell) : TAcceptability(G, community, shell) {
}

NetworKit::TDummyAcceptability::~TDummyAcceptability() {
}

double NetworKit::TDummyAcceptability::getValue(node v) {
	return 0.5;
}

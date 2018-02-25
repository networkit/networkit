/*
 * TopHarmonicCloseness.cpp
 *
 * Created on: 25.02.2018
 *		 Author: Eugenio Angriman
 */

#include "TopHarmonicCloseness.h"

namespace NetworKit {

TopHarmonicCloseness::TopHarmonicCloseness(const Graph &G, count k,
                                           bool first_heu, bool sec_heu)
    : G(G), k(k), first_heu(first_heu), sec_heu(sec_heu) {}

void TopHarmonicCloseness::run() {}
} // namespace NetworKit

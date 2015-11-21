/*
 * EliminationStage.cpp
 *
 *  Created on: 09.01.2015
 *      Author: Michael
 */

#include "EliminationStage.h"

namespace NetworKit {

EliminationStage::EliminationStage(const CSRMatrix &P, const Vector &q, const std::vector<index> &fSet, const std::vector<index> &cSet) : P(P), R(P.transpose()), q(q), fSet(fSet), cSet(cSet) {
}

const CSRMatrix& EliminationStage::getP() const {
	return P;
}

const CSRMatrix& EliminationStage::getR() const {
	return R;
}

const Vector& EliminationStage::getQ() const {
	return q;
}

const std::vector<index>& EliminationStage::getFSet() const {
	return fSet;
}

const std::vector<index>& EliminationStage::getCSet() const {
	return cSet;
}

count EliminationStage::getN() const {
	return fSet.size() + cSet.size();
}

} /* namespace NetworKit */

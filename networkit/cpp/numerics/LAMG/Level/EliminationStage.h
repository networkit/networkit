/*
 * EliminationStage.h
 *
 *  Created on: 09.01.2015
 *      Author: Michael
 */

#ifndef ELIMINATIONSTAGE_H_
#define ELIMINATIONSTAGE_H_

#include "../../../algebraic/CSRMatrix.h"
#include "../../../algebraic/Vector.h"

namespace NetworKit {

/**
 * @ingroup numerics
 */
class EliminationStage {
private:
	CSRMatrix P; // interpolation matrix
	CSRMatrix R;
	Vector q; // coarse result correction vector
	std::vector<index> fSet;
	std::vector<index> cSet;

public:
	EliminationStage(const CSRMatrix &P, const Vector &q, const std::vector<index> &fSet, const std::vector<index> &cSet);

	const CSRMatrix& getP() const;
	const CSRMatrix& getR() const;
 	const Vector& getQ() const;
	const std::vector<index>& getFSet() const;
	const std::vector<index>& getCSet() const;

	count getN() const;

};

} /* namespace NetworKit */

#endif /* ELIMINATIONSTAGE_H_ */

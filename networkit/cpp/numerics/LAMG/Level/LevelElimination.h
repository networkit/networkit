/*
 * LevelElimination.h
 *
 *  Created on: 10.01.2015
 *      Author: Michael
 */

#ifndef LEVELELIMINATION_H_
#define LEVELELIMINATION_H_

#include "Level.h"
#include "EliminationStage.h"

namespace NetworKit {

class LevelElimination : public Level {
private:
	std::vector<EliminationStage> coarseningStages;
	std::vector<index> cIndexFine;
	std::vector<Vector> bStages;

	void subVectorExtract(Vector &subVector, const Vector &vector, const std::vector<index> &elements) const;

public:
	LevelElimination(const CSRMatrix &A, const std::vector<EliminationStage> &coarseningStages);

	void coarseType(const Vector &xf, Vector &xc) const;
	void restrict(const Vector &bf, Vector &bc);
	void interpolate(const Vector &xc, Vector &xf) const;

	static count lvl;
};

} /* namespace NetworKit */

#endif /* LEVELELIMINATION_H_ */

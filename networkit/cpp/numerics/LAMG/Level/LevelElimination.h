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

/**
 * @ingroup numerics
 */
class LevelElimination : public Level {
private:
	std::vector<EliminationStage> coarseningStages;
	std::vector<index> cIndexFine;

	void subVectorExtract(Vector &subVector, const Vector &vector, const std::vector<index> &elements) const;

public:
	LevelElimination(const CSRMatrix &A, const std::vector<EliminationStage> &coarseningStages);

	void coarseType(const Vector &xf, Vector &xc) const;
	void restrict(const Vector &bf, Vector &bc, std::vector<Vector> &bStages) const;
	void interpolate(const Vector &xc, Vector &xf, const std::vector<Vector> &bStages) const;

	static count lvl;
};

} /* namespace NetworKit */

#endif /* LEVELELIMINATION_H_ */

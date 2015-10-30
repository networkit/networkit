/*
 * LevelHierarchy.h
 *
 *  Created on: 10.01.2015
 *      Author: Michael
 */

#ifndef LEVELHIERARCHY_H_
#define LEVELHIERARCHY_H_

#include "Level/Level.h"
#include "Level/LevelFinest.h"
#include "Level/LevelElimination.h"
#include "Level/LevelAggregation.h"

namespace NetworKit {

class LevelHierarchy {
private:
	std::vector<LevelType> levelType;
	std::vector<index> levelIndex;
	std::vector<LevelElimination> eliminationLevels;
	std::vector<LevelAggregation> aggregationLevels;
	LevelFinest finestLevel;


public:
	LevelHierarchy();

	void addFinestLevel(const CSRMatrix &A);
	void addEliminationLevel(const CSRMatrix &A, const std::vector<EliminationStage> &coarseningStages);
	void addAggregationLevel(const CSRMatrix &A, const CSRMatrix &P, const CSRMatrix &R);

	count size() const;
	LevelType getType(index levelIdx) const;
	Level& at(index levelIdx);
	double cycleIndex(index levelIdx);


};

} /* namespace NetworKit */

#endif /* LEVELHIERARCHY_H_ */

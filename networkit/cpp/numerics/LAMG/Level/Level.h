/*
 * Level.h
 *
 *  Created on: 09.01.2015
 *      Author: Michael
 */

#ifndef LEVEL_H_
#define LEVEL_H_

#include "../../../algebraic/CSRMatrix.h"

namespace NetworKit {

enum LevelType {FINEST, // original problem
	ELIMINATION, // lowdegree node elimination
	AGGREGATION, // aggregation of nodes with high affinity
	COARSEST // coarsest level
};

/**
 * @ingroup numerics
 * Abstract base class for an LAMG Level.
 */
class Level {
protected:
	LevelType type;
	CSRMatrix A;

public:
	Level(LevelType type);
	Level(LevelType type, const CSRMatrix &A);
	virtual ~Level() {}

	const CSRMatrix& getLaplacian() const;

	count getNumberOfNodes() const;

	virtual void coarseType(const Vector &xf, Vector &xc) const = 0;

	virtual void restrict(const Vector &bf, Vector &bc) const;

	virtual void restrict(const Vector &bf, Vector &bc, std::vector<Vector> &bStages) const;

	virtual void interpolate(const Vector &xc, Vector &xf) const;

	virtual void interpolate(const Vector &xc, Vector &xf, const std::vector<Vector> &bStages) const;


};

} /* namespace NetworKit */

#endif /* LEVEL_H_ */

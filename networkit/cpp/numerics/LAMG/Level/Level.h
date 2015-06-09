/*
 * Level.h
 *
 *  Created on: 09.01.2015
 *      Author: Michael
 */

#ifndef LEVEL_H_
#define LEVEL_H_

#include "../../../algebraic/Matrix.h"
#include "../../../algebraic/CSRMatrix.h"

namespace NetworKit {

enum LevelType {FINEST, // original problem
	ELIMINATION, // lowdegree node elimination
	AGGREGATION, // aggregation of nodes with high affinity
	COARSEST // coarsest level
};

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

	virtual void restrict(const Vector &bf, Vector &bc) = 0;

	virtual void interpolate(const Vector &xc, Vector &xf) const = 0;


};

} /* namespace NetworKit */

#endif /* LEVEL_H_ */

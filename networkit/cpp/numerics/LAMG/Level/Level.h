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
template<class Matrix>
class Level {
protected:
	LevelType type;
	Matrix A;

public:
	Level(LevelType type) : type(type) {}
	Level(LevelType type, const Matrix& A) : type(type), A(A) {}
	virtual ~Level() {}

	inline const Matrix& getLaplacian() const {
		return A;
	}

	inline count getNumberOfNodes() const {
		return A.numberOfRows();
	}

	virtual void coarseType(const Vector& /*xf*/, Vector& /*xc*/) const {}

	virtual void restrict(const Vector& /*bf*/, Vector& /*bc*/) const {}

	virtual void restrict(const Vector& /*bf*/, Vector& /*bc*/, std::vector<Vector>& /*bStages*/) const {}

	virtual void interpolate(const Vector& /*xc*/, Vector& /*xf*/) const {}

	virtual void interpolate(const Vector& /*xc*/, Vector& /*xf*/, const std::vector<Vector>& /*bStages*/) const {}
};



} /* namespace NetworKit */

#endif /* LEVEL_H_ */

/*
 * AlgebraicGlobals.h
 *
 *  Created on: May 31, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_CPP_ALGEBRAIC_ALGEBRAICGLOBALS_H_
#define NETWORKIT_CPP_ALGEBRAIC_ALGEBRAICGLOBALS_H_

#include "../Globals.h"

namespace NetworKit {

/** Represents a matrix entry s.t. matrix(row, column) = value */
struct Triplet {
	index row;
	index column;
	double value;
};

/** Floating point epsilon to use in comparisons. */
constexpr double FLOAT_EPSILON = 1e-9;

}



#endif /* NETWORKIT_CPP_ALGEBRAIC_ALGEBRAICGLOBALS_H_ */

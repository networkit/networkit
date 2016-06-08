/*
 * Semirings.h
 *
 *  Created on: May 31, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_CPP_ALGEBRAIC_SEMIRINGS_H_
#define NETWORKIT_CPP_ALGEBRAIC_SEMIRINGS_H_

#include <algorithm>

// *****************************************************
// 					Semiring Definitions
// *****************************************************

/**
 * @ingroup algebraic
 * codomain = (-infty, +infty)
 */
class ArithmeticSemiring {
public:
	ArithmeticSemiring() = default;
	virtual ~ArithmeticSemiring() = default;

	inline static double add(double a, double b) {
		return a + b;
	}

	inline static double mult(double a, double b) {
		return a * b;
	}

	inline static double zero() {return 0;};

	inline static double one() {return 1;};
};

/**
 * @ingroup algebraic
 * codomain = (-infty, +infty]
 */
class MinPlusSemiring {
public:
	MinPlusSemiring() = default;
	virtual ~MinPlusSemiring() = default;

	inline static double add(double a, double b) {
		return std::min(a,b);
	}

	inline static double mult(double a, double b) {
		return a+b;
	}

	inline static double zero() {return std::numeric_limits<double>::infinity();};

	inline static double one() {return 0;};
};

/**
 * @ingroup algebraic
 * codomain = [-infty, +infty)
 */
class MaxPlusSemiring {
public:
	MaxPlusSemiring() = default;
	virtual ~MaxPlusSemiring() = default;

	inline static double add(double a, double b) {
		return std::max(a,b);
	}

	inline static double mult(double a, double b) {
		return a+b;
	}

	inline static double zero() {return -std::numeric_limits<double>::infinity();};

	inline static double one() {return 0;};
};

/**
 * @ingroup algebraic
 * codomain = [0, +infty)
 */
class MinMaxSemiring {
public:
	MinMaxSemiring() = default;
	virtual ~MinMaxSemiring() = default;

	inline static double add(double a, double b) {
		return std::min(a,b);
	}

	inline static double mult(double a, double b) {
		return std::max(a,b);
	}

	inline static double zero() {return std::numeric_limits<double>::infinity();};

	inline static double one() {return 0;};
};

/**
 * @ingroup algebraic
 * codomain = (-infty, 0]
 */
class MaxMinSemiring{
public:
	MaxMinSemiring() = default;
	virtual ~MaxMinSemiring() = default;

	inline static double add(double a, double b) {
		return std::max(a,b);
	}

	inline static double mult(double a, double b) {
		return std::min(a,b);
	}

	inline static double zero() {return -std::numeric_limits<double>::infinity();};

	inline static double one() {return 0;};
};

/**
 * @ingroup algebraic
 */
class IntLogicalSemiring {
public:
	IntLogicalSemiring() = default;
	virtual ~IntLogicalSemiring() = default;

	inline static double add(double a, double b) {
		return (int) a || (int) b;
	}

	inline static double mult(double a, double b) {
		return (int) a && (int) b;
	}

	inline static double zero() {return 0;};

	inline static double one() {return 1;};
};

/**
 * @ingroup algebraic
 * codomain = [0, 1]
 */
class GaloisFieldSemiring {
public:
	GaloisFieldSemiring() = default;
	virtual ~GaloisFieldSemiring() = default;

	inline static double add(double a, double b) {
		return (int) a ^ (int) b;
	}

	inline static double mult(double a, double b) {
		return (int) a & (int) b;
	}

	inline static double zero() {return 0;};

	inline static double one() {return 1;};
};



#endif /* NETWORKIT_CPP_ALGEBRAIC_SEMIRINGS_H_ */

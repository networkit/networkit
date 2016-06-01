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

class MinPlusSemiring {
public:
	MinPlusSemiring() = default;
	virtual ~MinPlusSemiring() = default;

	inline static double add(double a, double b) {
		return std::min(a,b);
	}

	inline static double mult(double a, double b) {
		// check whether a or b are the zero element. If yes, return it to conform to the law a x zero() = zero().
		if (a == zero()) {
			return a;
		} else if (b == zero()) {
			return b;
		} else {
			return a+b;
		}
	}

	inline static double zero() {return std::numeric_limits<double>::max();};

	inline static double one() {return 0;};
};

class MaxPlusSemiring {
public:
	MaxPlusSemiring() = default;
	virtual ~MaxPlusSemiring() = default;

	inline static double add(double a, double b) {
		return std::max(a,b);
	}

	inline static double mult(double a, double b) {
		// check whether a or b are the zero element. If yes, return it to conform to the law a x zero() = zero().
		if (a == zero()) {
			return a;
		} else if (b == zero()) {
			return b;
		} else {
			return a+b;
		}
	}

	inline static double zero() {return -std::numeric_limits<double>::max();};

	inline static double one() {return 0;};
};

/**
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

	inline static double zero() {return std::numeric_limits<double>::max();};

	inline static double one() {return 0;};
};

/**
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

	inline static double zero() {return -std::numeric_limits<double>::max();};

	inline static double one() {return 0;};
};

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

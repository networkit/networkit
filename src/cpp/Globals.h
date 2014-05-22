/*
 * Globals.h
 *
 *  Created on: 06.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <cstdint>

namespace NetworKit {
	typedef uint64_t index; // more expressive name for an index into an array
	typedef uint64_t count; // more expressive name for an integer quantity
	typedef index node; // node indices are 0-based
}

extern bool PRINT_PROGRESS;
extern bool RAND_ORDER;
extern uint64_t INACTIVE_SEEDS;
extern bool NORMALIZE_VOTES;
extern double SCALE_STRENGTH;

extern uint64_t MIN_NUM_COMMUNITIES;
extern double REL_REPEAT_THRSH;

extern bool CALC_DISSIMILARITY;

extern int MAX_LOUVAIN_ITERATIONS;

#endif /* GLOBALS_H_ */

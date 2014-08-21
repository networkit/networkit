/*
 * Globals.h
 *
 *  Created on: 06.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <cstdint>
#include <cmath>
#include <limits>

namespace NetworKit {
	/** Typedefs **/

	typedef uint64_t index; // more expressive name for an index into an array
	typedef uint64_t count; // more expressive name for an integer quantity
	typedef index node; // node indices are 0-based
	typedef double edgeweight; // edge weight type
	typedef index edgeid;	// edge id

	/** Constants **/
	constexpr index none = std::numeric_limits<index>::max(); // value for not existing nodes/edges
	constexpr edgeweight defaultEdgeWeight = 1.0;
	constexpr edgeweight nullWeight = 0.0;
}

constexpr double PI = 2.0*std::acos(0);

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

/*
 * Parallelism.h
 *
 * Control functions related to parallelism.
 *
 *  Created on: 13.11.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef PARALLELISM_H_
#define PARALLELISM_H_



namespace Aux {

 // OpenMP
#ifdef _OPENMP
#include <omp.h>
inline void setNumberOfThreads(int nThreads) {
	omp_set_num_threads(nThreads);
}

inline int getCurrentNumberOfThreads() {
	return omp_get_num_threads();
}


inline int getMaxNumberOfThreads() {
	return omp_get_max_threads();
}

inline void enableNestedParallelism() {
	omp_set_nested(1); // enable nested parallelism
}

inline int getThreadNumber() {
	return omp_get_thread_num();
}

#else
// no OpenMP available

/**
 *  Set the number of threads available to the program.
 */
inline void setNumberOfThreads(int) {
}

/**
 * 
 * @return The number of threads currently running.
 */
inline int getCurrentNumberOfThreads() {
	return 1;
}

/**
 * 
 * @return The maximum number of threads available to the program.
 */
inline int getMaxNumberOfThreads() {
	return 1;
}

/** 
 * Enable OpenMP nested parallelism
 */
inline void enableNestedParallelism() {
}

/**
 * 
 * @return The current thread id.
 */
inline int getThreadNumber() {
	return 0;
}

#endif
}
#endif
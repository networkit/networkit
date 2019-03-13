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


/**
 *  Set the number of threads available to the program.
 */
void setNumberOfThreads(int nThreads);

/**
 *
 * @return The number of threads currently running.
 */
int getCurrentNumberOfThreads();


/**
 *
 * @return The maximum number of threads available to the program.
 */
int getMaxNumberOfThreads();


/** Enable OpenMP nested parallelism */
void enableNestedParallelism();

}

#endif

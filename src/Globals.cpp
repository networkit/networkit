/*
 * Globals.cpp
 *
 *  Created on: 06.02.2013
 *      Author: cls
 */

bool PRINT_PROGRESS = true;
bool RAND_ORDER = false; 	///< don't randomize LabelPropagation node order by default, since this seems to be faster and slightly better
bool NO_RECURSION = false; ///< don't recurse in ensemble clustering if set to true
bool NORMALIZE_VOTES = false; ///< if set, normalize votes of LP by dividing edge weight by weighted degree
bool SCALE_STRENGTH = false; ///< if set, scale cluster strengths in LP


int MIN_NUM_COMMUNITIES = 2; ///< minimum number of communities, so far only used by agglomerative clustering
double REL_REPEAT_THRSH = 5e-3; ///< threshold for minimum number of matching edges relative to number of vertices to proceed agglomeration

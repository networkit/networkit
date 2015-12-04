/*
 * LAMGSettings.h
 *
 *  Created on: 09.01.2015
 *      Author: Michael
 */

#ifndef LAMGSETTINGS_H_
#define LAMGSETTINGS_H_

namespace NetworKit {

// #initial test vectors (TVs) at each level.
#define TV_NUM 4
// #TVs to add upon each aggregation coarsening
#define TV_INC 1
// maximum allowed#TVs
#define TV_MAX 10
// #global sweeps to perform on each initial TV
#define SETUP_TV_SWEEPS 4
// max size for direct solver
#define MAX_DIRECT_SOLVE_SIZE 200
// maximum #aggregation coarsening levels to construct during the setup phase
#define SETUP_MAX_AGG_LEVELS 100
// maximum TOTAL # of coarsening levels to construct during the setup phase
#define SETUP_MAX_LEVELS 100
// solution cycle index. Also the design cycle index during setup phase.
#define SETUP_CYCLE_INDEX 1.5
// minimum number of sweeps to run to estimate relax ACF
#define SETUP_RELAX_ACF_MIN_SWEEPS 7
// if relaxation converges at this rate or faster, this becomes the CoarsestLevel
#define SETUP_MAX_COARSE_RELAX_ACF 0.3

#define MAX_COMBINED_ITERATES 4

/**************************
 * SETUP - Elimination    *
 **************************/

// maximum degree up to which a node gets eliminated
#define SETUP_ELIMINATION_MAX_DEGREE 4
// maximum number of elimination stages
#define SETUP_ELIMINATION_MAX_STAGES 1000
// node elimination stops if number of coarse nodes is only this fraction of the initial number of nodes
#define SETUP_ELIMINATION_MIN_ELIM_FRACTION 0.01

/**************************
 * SETUP - Aggregation    *
 **************************/

#define SETUP_AGGREGATION_WEAK_EDGE_THRESHOLD 0.1
// all nodes with a degree equal or higher than this threshold are marked as high-degree nodes
#define SETUP_AGGREGATION_DEGREE_THRESHOLD 8
// #sweeps (nu) to use during coarsening stages
#define SETUP_NU_DEFAULT 3

#define SETUP_COARSENING_WORK_GUARD 0.7

#define SETUP_MIN_AGGREGATION_STAGES 1

#define SETUP_MAX_AGGREGATION_STAGES 2

#define SETUP_RELAX_COARSEST_SWEEPS 400

}



#endif /* LAMGSETTINGS_H_ */

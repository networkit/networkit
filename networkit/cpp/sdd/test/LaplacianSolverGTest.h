/*
 * LaplacianSolverGTest.h
 *
 *  Created on: May 03, 2014
 *      Author: dhoske
 */

#ifndef NOGTEST
#ifndef LAPLACIAN_SOLVER_TEST_H_
#define LAPLACIAN_SOLVER_TEST_H_

#include <gtest/gtest.h>

#include "../LaplacianSolver.h"
#include "../SDDSolver.h"
#include "../SparsifySolver.h"

namespace NetworKit {
namespace SDD {

/** @internal */
class LaplacianSolverGTest : public testing::Test {
  /* EMPTY */
};
/** @endinternal */

}
}

#endif
#endif

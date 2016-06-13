/*
 * CSRMatrixGTest.h
 *
 *  Created on: May 13, 2015
 *      Author: Michael Wegner
 */

#ifndef CSRMATRIXGTEST_H_
#define CSRMATRIXGTEST_H_

#include "gtest/gtest.h"

#include "../CSRMatrix.h"
#include "../Vector.h"
#include "../../graph/Graph.h"
#include "../../io/METISGraphReader.h"

namespace NetworKit {

class CSRMatrixGTest : public testing::Test {
public:
	CSRMatrixGTest();
	virtual ~CSRMatrixGTest();
};

} /* namespace NetworKit */

#endif /* CSRMATRIXGTEST_H_ */

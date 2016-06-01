/*
 * GaussSeidelRelaxationGTest.h
 *
 *  Created on: 03.11.2014
 *      Author: Michael
 */

#ifndef NOGTEST

#ifndef GAUSSSEIDELRELAXATIONGTEST_H_
#define GAUSSSEIDELRELAXATIONGTEST_H_

#include "gtest/gtest.h"

#include "../../algebraic/CSRMatrix.h"
#include "../../algebraic/Vector.h"
#include "../GaussSeidelRelaxation.h"

namespace NetworKit {

class GaussSeidelRelaxationGTest : public testing::Test {
public:
	GaussSeidelRelaxationGTest() {}
	~GaussSeidelRelaxationGTest() {}
};

} /* namespace NetworKit */

#endif /* GAUSSSEIDELRELAXATIONGTEST_H_ */

#endif

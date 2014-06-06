/*
 * ConductanceTreeGTest.h
 *
 *  Created on: Jul 15, 2013
 *      Author: Henning
 */

#ifndef NOGTEST

#ifndef ConductanceTreeGTEST_H_
#define ConductanceTreeGTEST_H_

#include <gtest/gtest.h>

#include "../ConductanceTree.h"
#include "../../graph/Graph.h"
#include "../../graph/GraphGenerator.h"
#include "../../io/METISGraphReader.h"

namespace NetworKit {

class ConductanceTreeGTest: public testing::Test {
public:
	ConductanceTreeGTest();
	virtual ~ConductanceTreeGTest();
};

}

#endif /* ConductanceTreeGTEST_H_ */

#endif

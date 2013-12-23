#ifndef NOGTEST

#ifndef PRINTTESTHOSKE_H_
#define PRINTTESTHOSKE_H_

#include <gtest/gtest.h>

#include "../Graph.h"
#include "../../io/METISGraphReader.h"
#include "../../graph/Print_Hoske.h"

namespace NetworKit {

class PrintTest_Hoske: public testing::Test {
public:
	PrintTest_Hoske();
	virtual ~PrintTest_Hoske();
};

} /* namespace NetworKit */
#endif

#endif /*NOGTEST */

/*
 * SCDGTest.h
 *
 *  Created on: 06.06.2013
 *      Author: cls
 */

#ifndef NOGTEST


#ifndef SCDGTEST_H_
#define SCDGTEST_H_

#include <gtest/gtest.h>
#include "../../graph/Graph.h"
#include "../GreedyCommunityExpansion.h"

namespace NetworKit {

class SCDGTest: public testing::Test {
	public:
		SCDGTest();
		virtual ~SCDGTest();
};

} /* namespace NetworKit */
#endif /* SCDGTEST_H_ */

#endif /*NOGTEST */

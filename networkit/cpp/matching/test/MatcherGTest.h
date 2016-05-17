/*
 * MatcherGTest.h
 *
 *  Created on: Jun 14, 2013
 *      Author: Henning
 */

#ifndef NOGTEST

#ifndef MATCHERGTEST_H_
#define MATCHERGTEST_H_

#include <gtest/gtest.h>


namespace NetworKit {

class MatcherGTest: public testing::Test {
public:
	MatcherGTest() = default;
	virtual ~MatcherGTest() = default;
};

}

#endif /* MATCHERGTEST_H_ */

#endif

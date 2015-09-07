/*
 * GlobalGTest.h
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#ifndef NOGTEST

#ifndef GLOBALGTEST_H_
#define GLOBALGTEST_H_

#include <gtest/gtest.h>

namespace NetworKit {

class GlobalGTest: public testing::Test {
public:
	GlobalGTest();
	virtual ~GlobalGTest();
};

} /* namespace NetworKit */
#endif /* GLOBALGTEST_H_ */

#endif /* NOGTEST */

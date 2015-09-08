/*
 * GeneratorsGTest.h
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#ifndef NOGTEST

#ifndef GENERATORSGTEST_H_
#define GENERATORSGTEST_H_

#include <gtest/gtest.h>


namespace NetworKit {

class GeneratorsGTest: public testing::Test {
public:
	GeneratorsGTest();
	virtual ~GeneratorsGTest();
};

} /* namespace NetworKit */
#endif /* GENERATORSGTEST_H_ */

#endif /*NOGTEST */

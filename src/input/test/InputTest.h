/*
 * InputTest.h
 *
 *  Created on: 06.12.2012
 *      Author: cls
 */

#ifndef InputTest_H_
#define InputTest_H_

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "../../aux/log.h"
#include "../METISParser.h"

namespace EnsembleClustering {

class InputTest : public CPPUNIT_NS::TestFixture {

	CPPUNIT_TEST_SUITE(InputTest);
	CPPUNIT_TEST(testMETISParser);
	CPPUNIT_TEST_SUITE_END();

protected:

	std::string graphFilePath;

public:

	InputTest();

	virtual ~InputTest();

	void setUp();

	void tearDown();

	/*** Tests **/

	void testMETISParser();

};

} /* namespace EnsembleClustering */
#endif /* InputTest_H_ */

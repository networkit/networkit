/*
 * InputTest.cpp
 *
 *  Created on: 06.12.2012
 *      Author: cls
 */

#include "InputTest.h"

namespace EnsembleClustering {

CPPUNIT_TEST_SUITE_REGISTRATION(InputTest);

InputTest::InputTest() {
	// TODO Auto-generated constructor stub

}

InputTest::~InputTest() {
	// TODO Auto-generated destructor stub
}

void InputTest::setUp() {
	this->graphFilePath = "/Users/cls/workspace/Data/DIMACS/kron_g500-simple-logn16.graph";
}

void InputTest::tearDown() {
}

void InputTest::testMETISParser() {
	bool success = false;
	METISParser* parser = new METISParser();

	parser->open(this->graphFilePath);
	std::pair<int, int> header = parser->getHeader();

	int lc = 0;
	while (parser->hasNext()) {
		parser->getNext();
		lc++;
	}

	LOG4CXX_INFO(log4cxx::Logger::getRootLogger(), "parsed " << lc << " lines");

	parser->close();
	success = true;
	CPPUNIT_ASSERT_EQUAL(success, true);

}

} /* namespace EnsembleClustering */

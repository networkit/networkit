/*
 * InputGTest.h
 *
 *  Created on: 12.12.2012
 *      Author: cls
 */

#ifndef INPUTGTEST_H_
#define INPUTGTEST_H_

#include <gtest/gtest.h>
#include <fstream>

#include "../../aux/Log.h"
#include "../METISParser.h"

#include "../GraphIO.h"
#include "../../graph/GraphGenerator.h"

namespace EnsembleClustering {

class InputGTest: public testing::Test {

};



//TEST_F(InputGTest, testMETISParser) {
//	bool success = false;
//	METISParser parser;
//
//	std::string graphFilePath = "/Users/cls/workspace/Data/DIMACS/kron_g500-simple-logn16.graph";
//
//	parser.open(graphFilePath);
//	std::pair<int, int> header = parser.getHeader();
//
//	int lc = 0;
//	while (parser.hasNext()) {
//		parser.getNext();
//		lc++;
//	}
//
//	parser.close();
//	success = true;
//	EXPECT_EQ(success, true);
//}

} /* namespace EnsembleClustering */
#endif /* INPUTGTEST_H_ */

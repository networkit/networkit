/*
 * DynamicsGTest.cpp
 *
 *  Created on: 24.12.2013
 *      Author: cls
 */

#include "DynamicsGTest.h"

#include "../DGSStreamParser.h"

namespace NetworKit {

TEST_F(DynamicsGTest, testDGSStreamParser) {
	DGSStreamParser parser("input/example2.dgs");
	auto stream = parser.getStream();
	for (auto event : stream) {
		std::cout << event.toString() << " ";
	}
	std::cout << std::endl;
}

} /* namespace NetworKit */

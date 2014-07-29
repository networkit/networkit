/*
 * ApproxBetweennessGTest.cpp
 *
 *  Created on: 30.06.2014
 *      Author: moritzl
 */

#include "ApproxBetweennessGTest.h"
#include "../ApproxBetweenness.h"
#include "../../generators/ErdosRenyiGenerator.h"

namespace NetworKit {

ApproxBetweennessGTest::ApproxBetweennessGTest() {
	// TODO Auto-generated constructor stub

}

ApproxBetweennessGTest::~ApproxBetweennessGTest() {
	// TODO Auto-generated destructor stub
}

TEST_F(ApproxBetweennessGTest, testApproxDiameterErdos) {
	ErdosRenyiGenerator gen(10000,0.001);
	Graph G1 = gen.generate();
	ApproxBetweenness approx(G1, 0.05, 0.1, 20);
	approx.run();
}

}

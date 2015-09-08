/*
 * GlobalGTest.cpp
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#ifndef NOGTEST

#include "GlobalGTest.h"

#include "../ClusteringCoefficient.h"

#include "../../generators/ErdosRenyiGenerator.h"

namespace NetworKit {

GlobalGTest::GlobalGTest() {

}

GlobalGTest::~GlobalGTest() {

}


TEST_F(GlobalGTest, testClusteringCoefficient) {

	ErdosRenyiGenerator graphGen(10, 1.0);
	Graph G = graphGen.generate();

	ClusteringCoefficient clusteringCoefficient;
	double cc = clusteringCoefficient.avgLocal(G);

	EXPECT_EQ(1.0, cc);
}









} /* namespace NetworKit */

#endif /*NOGTEST*/

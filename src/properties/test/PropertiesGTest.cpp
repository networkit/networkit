/*
 * PropertiesGTest.cpp
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#include "PropertiesGTest.h"

namespace NetworKit {

PropertiesGTest::PropertiesGTest() {
	// TODO Auto-generated constructor stub

}

PropertiesGTest::~PropertiesGTest() {
	// TODO Auto-generated destructor stub
}


TEST_F(PropertiesGTest, testClusteringCoefficient) {

	GraphGenerator gen;
	Graph G = gen.makeErdosRenyiGraph(100, 1.0);

	ClusteringCoefficient clusteringCoefficient;
	double cc = clusteringCoefficient.calculate(G);

	EXPECT_EQ(1.0, cc);

}


} /* namespace NetworKit */

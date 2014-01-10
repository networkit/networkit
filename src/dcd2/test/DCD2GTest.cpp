/*
 * DCD2GTest.cpp
 *
 *  Created on: 07.01.2014
 *      Author: cls
 */

#include "DCD2GTest.h"

#include "../DynPLP.h"
#include "../DynPLM.h"
#include "../../dynamics/DGSStreamParser.h"
#include "../../dynamics/GraphUpdater.h"

namespace NetworKit {


TEST_F(DCD2GTest, testGraphUpdater) {
	DGSStreamParser parser("input/example2.dgs", true, 1);
	DEBUG("getting event stream");
	std::vector<GraphEvent> stream = parser.getStream();
	Graph G;
	GraphUpdater gu(G);
	DEBUG("updating graph");
	gu.update(stream);

	EXPECT_EQ(3, G.numberOfNodes());

}


TEST_F(DCD2GTest, testDynPLP) {
	DGSStreamParser parser("input/example2.dgs");
	std::vector<GraphEvent> stream = parser.getStream();
	Graph G;
	DynPLP dynPLP(G, 0);

	GraphUpdater gu(G);
	gu.update(stream);

	dynPLP.process(stream);
	Clustering zeta = dynPLP.retrieve();

	EXPECT_TRUE(zeta.isProper(G));

}


TEST_F(DCD2GTest, testDynPLM) {
	DGSStreamParser parser("input/example2.dgs");
	std::vector<GraphEvent> stream = parser.getStream();
	Graph G;
	DynPLM dynPLM(G, 0);

	GraphUpdater gu(G);
	gu.update(stream);

	dynPLM.process(stream);
	Clustering zeta = dynPLM.retrieve();

	EXPECT_TRUE(zeta.isProper(G));
}





} /* namespace NetworKit */

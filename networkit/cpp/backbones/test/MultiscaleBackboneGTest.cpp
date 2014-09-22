/*
 * MultiscaleBackboneGTest.cpp
 *
 *  Created on: 20.06.2014
 *      Author: Gerd Lindner
 */

#ifndef NOGTEST

#include "MultiscaleBackboneGTest.h"

#include "../Backbones.h"
#include "../MultiscaleAttributizer.h"


namespace NetworKit {

TEST_F(MultiscaleBackboneGTest, testSimpleMultiscaleBackbone) {
	Graph g(6, true, false);

	g.setWeight(0, 1, 1.0);
	g.setWeight(0, 2, 5.0);
	g.setWeight(0, 3, 10);
	g.setWeight(0, 4, 20.0);
	g.setWeight(4, 5, 1.0);
	g.setWeight(3, 5, 0.5);
	g.indexEdges();

	//TODO: weight --> attribute

	MultiscaleAttributizer attributizer;
	EXPECT_NEAR(0.0878244, attributizer.getProbability(4, 0.5555), 1e-5) << "faulty probability calculation";
	/**
	 * a01 = 0.91896
	 * a02 = 0.639212
	 * a03 = 0.376716
	 * a04 = 0.0878244
	 * a45 = 0.952381
	 * a40 = 0.047619
	 * a35 = 0.952381
	 * a30 = 0.047619
	 * a54 = 0.33333333
	 * a53 = 0.66666666
	 */

	//Compare the backbone graph to the expected backbone.

	MultiscaleBackbone backbone(0.5);
	Graph b = backbone.calculate(g);

	EXPECT_EQ(3, b.numberOfEdges());
	EXPECT_TRUE(b.hasEdge(0, 4));
	EXPECT_TRUE(b.hasEdge(0, 3));
	EXPECT_TRUE(b.hasEdge(4, 5));

	MultiscaleBackbone backbone2(0.3333);
	b = backbone2.calculate(g);
	EXPECT_EQ(2, b.numberOfEdges());
	EXPECT_TRUE(b.hasEdge(0, 4));
	EXPECT_TRUE(b.hasEdge(0, 3));
}

}
/* namespace NetworKit */

#endif /*NOGTEST */

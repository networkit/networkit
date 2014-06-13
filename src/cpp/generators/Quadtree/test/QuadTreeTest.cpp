/*
 * QuadTreeTest.cpp
 *
 *  Created on: 28.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#include "QuadTreeTest.h"
#include "../../../auxiliary/Random.h"
#include "../../../auxiliary/Log.h"
#include "../../HyperbolicSpace.h"

namespace NetworKit {

QuadTreeTest::QuadTreeTest() {
	// TODO Auto-generated constructor stub

}

QuadTreeTest::~QuadTreeTest() {
	// TODO Auto-generated destructor stub
}

TEST_F(QuadTreeTest, testQuadTreeInsertion) {
	count n = 10000;
	double R = acosh((double)n/(2*M_PI)+1);
	vector<double> angles(n);
	vector<double> radii(n);
	HyperbolicSpace::fillPoints(&angles, &radii, 1, 1);
	Quadtree<index> quad(R);

	for (index i = 0; i < n; i++) {
		EXPECT_GE(angles[i], 0);
		EXPECT_LT(angles[i], 2*M_PI);
		EXPECT_GE(radii[i], 0);
		EXPECT_LT(radii[i], R);
		TRACE("Added (", angles[i], ",", radii[i], ")");
		quad.addContent(i, angles[i], radii[i]);
	}
	vector<index> all = quad.getElements();
	EXPECT_EQ(all.size(), n);
	index comparison = 0;
	DEBUG("Using ", comparison, " at (", angles[comparison], ",", radii[comparison], ") as query point");
	vector<index> closeToOne = quad.getCloseElements(angles[comparison], radii[comparison], R);
	EXPECT_LE(closeToOne.size(), n);

	for (index i = 0; i < closeToOne.size(); i++) {
		//no corrupt indices
		ASSERT_LE(closeToOne[i], n);

		//close results should actually be close
		EXPECT_LE(HyperbolicSpace::getDistance(angles[comparison], radii[comparison], angles[closeToOne[i]], radii[closeToOne[i]]), R);
		for (index j = 0; j < i; j++) {
			/**
			 * results are unique
			 */
			EXPECT_NE(closeToOne[i], closeToOne[j]);
		}
	}

	for (index i = 0; i < n; i++) {
		if (HyperbolicSpace::getDistance(angles[comparison], radii[comparison], angles[i], radii[i]) < R) {
			bool found = false;
			QuadNode<index> responsibleNode = * getRoot(quad).getAppropriateLeaf(angles[i], radii[i]);

			TRACE("Getting lower bound for responsible node");
			double bound = responsibleNode.distanceLowerBound(angles[comparison], radii[comparison]);
			double actualDistance = HyperbolicSpace::getDistance(angles[comparison], radii[comparison], angles[i], radii[i]);
			EXPECT_GE(actualDistance, bound);
			EXPECT_TRUE(responsibleNode.responsible(angles[i], radii[i]));

			for (index j = 0; j < closeToOne.size(); j++) {
				if (closeToOne[j] == i) {
					found = true;
					break;
				}
			}
			EXPECT_TRUE(found) << "dist(" << i << "," << comparison << ") = "
					<< HyperbolicSpace::getDistance(angles[comparison], radii[comparison], angles[i], radii[i]) << " < " << R;
			if (!found || actualDistance < bound) {
				DEBUG("angle: ", angles[i], ", radius: ", radii[i], ", leftAngle: ", responsibleNode.getLeftAngle(),
						", rightAngle: ", responsibleNode.getRightAngle(), ", minR: ", responsibleNode.getMinR(), ", maxR:", responsibleNode.getMaxR());
				DEBUG("dist(", comparison, ", leftMin)=", HyperbolicSpace::getDistance(angles[comparison], radii[comparison], responsibleNode.getLeftAngle(), responsibleNode.getMinR()));
				DEBUG("dist(", comparison, ", leftMax)=", HyperbolicSpace::getDistance(angles[comparison], radii[comparison], responsibleNode.getLeftAngle(), responsibleNode.getMaxR()));
				DEBUG("dist(", comparison, ", rightMin)=", HyperbolicSpace::getDistance(angles[comparison], radii[comparison], responsibleNode.getRightAngle(), responsibleNode.getMinR()));
				DEBUG("dist(", comparison, ", rightMax)=", HyperbolicSpace::getDistance(angles[comparison], radii[comparison], responsibleNode.getRightAngle(), responsibleNode.getMaxR()));
			}
		}
	}
}

} /* namespace NetworKit */

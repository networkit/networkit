/*
 * QuadTreeTest.cpp
 *
 *  Created on: 28.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#include "QuadTreeTest.h"
#include "../Quadtree.h"
#include "../../../auxiliary/Random.h"
#include "../../HyperbolicSpace.h"

namespace NetworKit {

QuadTreeTest::QuadTreeTest() {
	// TODO Auto-generated constructor stub

}

QuadTreeTest::~QuadTreeTest() {
	// TODO Auto-generated destructor stub
}

TEST_F(QuadTreeTest, testQuadTreeInsertion) {
	count n = 200;
	double R = acosh((double)n/(2*M_PI)+1);
	vector<double> angles(n);
	vector<double> radii(n);
	HyperbolicSpace::fillPoints(&angles, &radii, 1, 1);
	Quadtree<index> quad(R);

	for (index i = 0; i < n; i++) {
		quad.addContent(i, angles[i], radii[i]);
	}
	vector<index> all = quad.getElements();
	EXPECT_EQ(all.size(), n);
	index comparison = 0;
	vector<index> closeToOne = quad.getCloseElements(angles[comparison], radii[comparison], R);
	EXPECT_LE(closeToOne.size(), n);

	for (index i = 0; i < closeToOne.size(); i++) {
		ASSERT_LE(closeToOne[i], n);
		EXPECT_LE(HyperbolicSpace::getDistance(angles[comparison], radii[comparison], angles[closeToOne[i]], radii[closeToOne[i]]), R);
		for (index j = 0; j < i; j++) {
			EXPECT_NE(closeToOne[i], closeToOne[j]);
		}
	}

	for (index i = 0; i < n; i++) {
		if (HyperbolicSpace::getDistance(angles[comparison], radii[comparison], angles[i], radii[i]) < R) {
			bool found = false;
			for (index j = 0; j < closeToOne.size(); j++) {
				if (closeToOne[j] == i) {
					found = true;
					break;
				}
			}
			EXPECT_TRUE(found) << "dist(" << i << "," << comparison << ") = "
					<< HyperbolicSpace::getDistance(angles[comparison], radii[comparison], angles[i], radii[i]) << " < " << R;
			//<< "Node " << i << " at (" << angles[i] << "," << radii[i] << ") is close to node "
			//		<< comparison << " at (" << angles[comparison] << "," << radii[comparison] << "), but doesn't show up in list of size " << closeToOne.size();
		}
	}
}

} /* namespace NetworKit */

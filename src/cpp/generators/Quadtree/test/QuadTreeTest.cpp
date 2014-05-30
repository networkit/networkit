/*
 * QuadTreeTest.cpp
 *
 *  Created on: 28.05.2014
 *      Author: moritz
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
	Quadtree<index> quad(R);

	for (index i = 0; i < n; i++) {
		angles[i] = Aux::Random::real(0, 2*M_PI);
		/**
		 * for the radial coordinate distribution, I took the probability density from Greedy Forwarding in Dynamic Scale-Free Networks Embedded in Hyperbolic Metric Spaces
		 * f (r) = sinh r/(cosh R − 1)
		 * \int sinh = cosh+const
		 */

		double maxcdf = cosh(R);
		double random = Aux::Random::real(1, maxcdf);
		radii[i] = acosh(random);
		assert(radii[i] <= R);
		quad.addContent(i, angles[i], radii[i]);
	}
	vector<index> all = quad.getElements();
	EXPECT_EQ(all.size(), n);
	index comparison = 0;
	vector<index> closeToOne = quad.getCloseElements(angles[comparison], radii[comparison], R);
	EXPECT_LE(closeToOne.size(), n);

	HyperbolicSpace space(R);

	for (index i = 0; i < closeToOne.size(); i++) {
		EXPECT_LE(space.getDistance(angles[comparison], radii[comparison], angles[closeToOne[i]], radii[closeToOne[i]]), R);
	}
}

} /* namespace NetworKit */

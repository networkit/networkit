/*
 * QuadTreePolarEuclidGTest.cpp
 *
 *  Created on: 22.07.2015
 *      Author: moritzl
 */

#include <vector>
#include <cmath>
#include "../../../auxiliary/Random.h"
#include "QuadTreePolarEuclidGTest.h"
#include "../QuadtreePolarEuclid.h"

using std::vector;

namespace NetworKit {

QuadTreePolarEuclidGTest::QuadTreePolarEuclidGTest() {
	// TODO Auto-generated constructor stub

}

QuadTreePolarEuclidGTest::~QuadTreePolarEuclidGTest() {
	// TODO Auto-generated destructor stub
}

TEST_F(QuadTreePolarEuclidGTest, testQuadTreePolarEuclidInsertion) {
	/**
	 * setup of data structures and constants
	 */
	double maxR = 2;
	count n = 1000;
	vector<double> angles(n);
	vector<double> radii(n);
	vector<index> content(n);

	double minPhi = 0;
	double maxPhi = 2*M_PI;
	double minR = 0;

	/**
	 * get random number generators
	 */

	std::uniform_real_distribution<double> phidist{minPhi, maxPhi};
	std::uniform_real_distribution<double> rdist{minR, maxR};

	/**
	 * fill vectors
	 */
	for (index i = 0; i < n; i++) {
		angles[i] = phidist(Aux::Random::getURNG());
		radii[i] = rdist(Aux::Random::getURNG());
		content[i] = i;
	}

	QuadtreePolarEuclid<index> tree(angles, radii, content);
	EXPECT_EQ(n, tree.size());

	/**
	 * elements are returned
	 */
	vector<index> returned = tree.getElements();
	EXPECT_EQ(n, returned.size());
	sort(returned.begin(), returned.end());
	for (index i = 0; i < returned.size(); i++) {
		EXPECT_EQ(i, returned[i]);
	}
}


} /* namespace NetworKit */

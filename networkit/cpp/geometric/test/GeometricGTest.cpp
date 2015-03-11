/*
 * GeometricGTest.cpp
 *
 *  Created on: 29.07.2014
 *      Author: moritzl
 */

#include "GeometricGTest.h"

namespace NetworKit {

GeometricGTest::GeometricGTest() {
	// TODO Auto-generated constructor stub

}

GeometricGTest::~GeometricGTest() {
	// TODO Auto-generated destructor stub
}

TEST_F(GeometricGTest, testConversion) {
	Point2D<double> a(6,7);
	double angle, radius;
	HyperbolicSpace::cartesianToPolar(a, angle, radius);
	Point2D<double> back = HyperbolicSpace::polarToCartesian(angle, radius);
	EXPECT_LE(a[0] - back[0], 0.000001);
	EXPECT_LE(a[1] - back[1], 0.000001);
	count n = 1000;
	vector<double> angles(n);
	vector<double> radii(n);
	HyperbolicSpace::fillPoints(angles, radii, 1, 1);
	for (index i = 0; i < n; i++) {
		Point2D<double> point = HyperbolicSpace::polarToCartesian(angles[i], radii[i]);
		double phi, r;
		HyperbolicSpace::cartesianToPolar(point, phi,r);
		EXPECT_GE(phi, 0) << "Point (" << point[0] << "," << point[1] << ") was not converted correctly";
		EXPECT_GE(r, 0);
		EXPECT_LE(abs(phi - angles[i]), 0.000001);
		EXPECT_LE(abs(r - radii[i]), 0.000001);
	}
}

} /* namespace NetworKit */

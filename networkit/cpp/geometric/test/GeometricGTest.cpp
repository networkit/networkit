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


TEST_F(GeometricGTest, testIntersect) {
	Point2D<double> a(1,0);
	Point2D<double> b(1,1);
	Point2D<double> c(5,0);
	Point2D<double> d(0,2);
	Point2D<double> intersect = HyperbolicSpace::intersect(a, b, c, d);
	EXPECT_EQ(5, intersect[0]);
	EXPECT_EQ(4, intersect[1]);

}

TEST_F(GeometricGTest, testIntersectCircle) {
	Point2D<double> a(-5,0);
	Point2D<double> b(5,0);
	Point2D<double> c(0,5);
	Point2D<double> mid1 = (a+b).scale(0.5);
	Point2D<double> mid2 = (b+c).scale(0.5);
	EXPECT_EQ(0, mid1[0]);
	EXPECT_EQ(0, mid1[1]);
	EXPECT_EQ(2.5, mid2[0]);
	EXPECT_EQ(2.5, mid2[1]);

	Point2D<double> om1 = HyperbolicSpace::orth(a-b);
	Point2D<double> om2 = HyperbolicSpace::orth(b-c);
	EXPECT_EQ(0, om1[0]);
	EXPECT_EQ(-10, om1[1]);
	EXPECT_EQ(5, om2[0]);
	EXPECT_EQ(5, om2[1]);
	Point2D<double> intersect = HyperbolicSpace::intersect(mid1, om1, mid2, om2);

	EXPECT_EQ(0, intersect[0]);
	EXPECT_EQ(0, intersect[1]);

	//TODO: test other cases
}

TEST_F(GeometricGTest, testMirrorOnCircle) {
	Point2D<double> a(5,0);
	Point2D<double> origin(0,0);
	double radius = 6;
	Point2D<double> image = HyperbolicSpace::mirrorOnCircle(a, origin, radius);

	EXPECT_LE(radius*radius/5-image[0], 0.000001);
	EXPECT_EQ(0, image[1]);
}

TEST_F(GeometricGTest, testCoordinateTransformation) {
	Point2D<double> a(5,2);
	Point2D<double> b(4,1);
	Point2D<double> m;
	double radius;
	double bound = 10;
	HyperbolicSpace::getTransmutationCircle(a, b, bound, m, radius);
	DEBUG("Circle center at (",m[0], ",", m[1], ") with radius ", radius);
	DEBUG("d(a,m)=", a.distance(m));
	DEBUG("d(b,m)=", b.distance(m));
	EXPECT_GE(m.length(), bound);
	EXPECT_NE(m.distance(a), radius); //for mirroring, none of the points can be on the circle
	EXPECT_NE(m.distance(b), radius);
	if (m.distance(a) < radius) EXPECT_GE(m.distance(b), radius); //the circle has to be between the points
	else if (m.distance(a) > radius) EXPECT_LE(m.distance(b), radius);

	Point2D<double> mirrored = HyperbolicSpace::mirrorOnCircle(a, m, radius);
	DEBUG("Mirrored a(", a[0], ",", a[1], ") to (", mirrored[0], ",", mirrored[1], ")");
	DEBUG("d(mirror,m)=", mirrored.distance(m));
	EXPECT_LE(mirrored.distance(b), 0.00001);
}

TEST_F(GeometricGTest, testCircleCenter) {
	Point2D<double> a(-5,0);
	Point2D<double> b(5,0);
	Point2D<double> c(0,5);
	Point2D<double> center = HyperbolicSpace::circleCenter(a, b, c);
	EXPECT_EQ(0, center[0]);
	EXPECT_EQ(0, center[1]);
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
	HyperbolicSpace::fillPoints(&angles, &radii, 1, 1);
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



TEST_F(GeometricGTest, testIsometries) {
	Point2D<double> a(0.5,0);
	Point2D<double> b(0.7,0);
	double dist =  acosh(  1 + 2*a.squaredDistance(b) / ((1 - a.squaredLength())*(1 - b.squaredLength()))  );
	//EXPECT_EQ(0.2, HyperbolicSpace::getHyperbolicDistance(a,b));

	Point2D<double> c(0.4,0);
	Point2D<double> origin(0,0);
	double R = 1;
	Point2D<double> m;
	double radius;

	HyperbolicSpace::getTransmutationCircle(c, origin, R, m, radius);
	DEBUG("Circle center at (",m[0], ",", m[1], ") with radius ", radius);
	EXPECT_LE(radius, m.length());
	EXPECT_GE(radius, m.distance(c));
	EXPECT_EQ(radius*radius+R*R, m.length()*m.length());
	Point2D<double> adash = HyperbolicSpace::mirrorOnCircle(a, m, radius);
	Point2D<double> bdash = HyperbolicSpace::mirrorOnCircle(b, m, radius);
	EXPECT_LE(adash.length() , R);
	EXPECT_LE(bdash.length() , R);
	DEBUG("Mirrored a(", a[0], ",", a[1], ") to (", adash[0], ",", adash[1], ")");
	DEBUG("Mirrored b(", b[0], ",", b[1], ") to (", bdash[0], ",", bdash[1], ")");
	double distdashnom = 2*adash.squaredDistance(bdash);
	double distdashdenom = (1 - adash.squaredLength())*(1 - bdash.squaredLength());
	double distdash = acosh(  1 +  distdashnom/  distdashdenom );
	EXPECT_LE(dist - distdash, 0.0001);
}

TEST_F(GeometricGTest, testPointOnCircle) {
	count n = 1000;
	Point2D<double> origin;
	vector<double> angles(n);
	vector<double> radii(n);
	HyperbolicSpace::fillPoints(&angles, &radii, 1, 1);
	for (index i = 0; i < n; i++) {
		Point2D<double> a = HyperbolicSpace::polarToCartesian(angles[i], radii[i]);
		double radius = Aux::Random::real(HyperbolicSpace::getHyperbolicDistance(origin, a));
		Point2D<double> second = HyperbolicSpace::getPointOnHyperbolicCircle(a, radius);
		EXPECT_LE(abs(radius-HyperbolicSpace::getHyperbolicDistance(a, second)), 0.0001);
	}
}

TEST_F(GeometricGTest, testEuclideanCircleProjection) {
	count n = 1000;
	Point2D<double> origin;
		vector<double> angles(n);
		vector<double> radii(n);
		HyperbolicSpace::fillPoints(&angles, &radii, 1, 1);
		for (index i = 0; i < n; i++) {
			/**
			 * get two hyperbolic points
			 */
			Point2D<double> a = HyperbolicSpace::polarToCartesian(angles[i], radii[i]);
			double radius = Aux::Random::real(HyperbolicSpace::getHyperbolicDistance(origin, a));
			Point2D<double> second = HyperbolicSpace::getPointOnHyperbolicCircle(a, radius);

			//get center and radius of euclidean circle
			Point2D<double> center;
			double euRadius;
			HyperbolicSpace::getEuclideanCircle(a, second, center, euRadius);
			EXPECT_LE(euRadius + center.length(), 1);

			//get highest point and check
			Point2D<double> highest = center;
			highest.scale((center.length()+euRadius)/center.length());
			EXPECT_LE(abs(radius-HyperbolicSpace::getHyperbolicDistance(a, highest)), 0.0001);

			//get lowest point and check
			Point2D<double> lowest = center;
			lowest.scale((center.length()-euRadius)/center.length());
			EXPECT_LE(abs(radius-HyperbolicSpace::getHyperbolicDistance(a, lowest)), 0.0001);

			//check whether center lies on arc between a and origin
			EXPECT_LE(HyperbolicSpace::getHyperbolicDistance((highest+lowest).scale(0.5),a), radius);
			double phi_a, r_a, phi_c, r_c;
			HyperbolicSpace::cartesianToPolar(a, phi_a, r_a);
			HyperbolicSpace::cartesianToPolar(center, phi_c, r_c);

			EXPECT_LE(abs(phi_c - phi_a), 0.000001);

			//numeric solution:
			double r_c2 = (2*r_a)/((1-r_a*r_a)*(cosh(radius)-1+2/(1-r_a*r_a)));
			EXPECT_LE(abs(r_c-r_c2), 0.000001);

			double euRadius2 = sqrt(r_c2*r_c2 - (2*r_a*r_a - (1-r_a*r_a)*(cosh(radius)-1))/((1-r_a*r_a)*(cosh(radius)-1+2/(1-r_a*r_a))));
			EXPECT_LE(abs(euRadius-euRadius2), 0.000001);
		}
}

} /* namespace NetworKit */

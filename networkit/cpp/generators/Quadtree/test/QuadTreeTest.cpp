/*
 * QuadTreeTest.cpp
 *
 *  Created on: 28.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#include "QuadTreeTest.h"
#include "../../../auxiliary/Random.h"
#include "../../../auxiliary/Log.h"
#include "../../../geometric/HyperbolicSpace.h"

namespace NetworKit {

QuadTreeTest::QuadTreeTest() {
	// TODO Auto-generated constructor stub

}

QuadTreeTest::~QuadTreeTest() {
	// TODO Auto-generated destructor stub
}

TEST_F(QuadTreeTest, testQuadTreeInsertion) {
	count n = 1000;
	double R = 1;
	vector<double> angles(n);
	vector<double> radii(n);
	HyperbolicSpace::fillPoints(&angles, &radii, 1, 1);
	double max = 0;
	for (index i = 0; i < n; i++) {
		if (radii[i] > max) {
			max = radii[i];
		}
	}
	Quadtree<index> quad(max+(1-max)/4);

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
	for (index testindex = 0; testindex < 100; testindex++) {
		index comparison = Aux::Random::integer(n-1);
		Point2D<double> origin;
		Point2D<double> query = HyperbolicSpace::polarToCartesian(angles[comparison], radii[comparison]);
		DEBUG("Using ", comparison, " at (", angles[comparison], ",", radii[comparison], ") as query point");

		vector<index> closeToOne = quad.getCloseElements(HyperbolicSpace::polarToCartesian(angles[comparison], radii[comparison]), R);
		EXPECT_LE(closeToOne.size(), n);

		/**
		* probable query circle

		Point<double> pointOnEdge = HyperbolicSpace::getPointOnHyperbolicCircle(query, R);
		double distance = HyperbolicSpace::getHyperbolicDistance(query, pointOnEdge);
		EXPECT_LE(abs(distance - R), 0.00001);
		Point<double> center;
		double radius, minPhi, maxPhi;
		HyperbolicSpace::getEuclideanCircle(query, pointOnEdge, center, radius);
		DEBUG("Assuming circle at (", center[0], ",",center[1], ") with radius ", radius);
		double minR = center.length() - radius;
		double maxR = center.length() + radius;
		assert(maxR < 1);
		if (minR < 0) {
			maxR = std::max(abs(minR), maxR);
			minR = 0;
			minPhi = 0;
			maxPhi = 2*M_PI;
		} else {
			double spread = asin(radius / center.length());
			double phi_c, r_c;
			HyperbolicSpace::cartesianToPolar(center, phi_c, r_c);
			minPhi = phi_c - spread;
			maxPhi = phi_c + spread;

			 //what to do if they overlap the 2\pi line? Well, have to make two separate calls and collect
			/
		}
		*/
		for (index i = 0; i < closeToOne.size(); i++) {
			//no corrupt indices
			ASSERT_LE(closeToOne[i], n);

			//close results should actually be close
			EXPECT_LE(HyperbolicSpace::getHyperbolicDistance(angles[comparison], radii[comparison], angles[closeToOne[i]], radii[closeToOne[i]]), R);
			for (index j = 0; j < i; j++) {
				/**
				 * results are unique
				 */
				EXPECT_NE(closeToOne[i], closeToOne[j]);
			}
		}
		count notfound = 0;
		count didfind = 0;
		for (index i = 0; i < n; i++) {
			if (HyperbolicSpace::getHyperbolicDistance(angles[comparison], radii[comparison], angles[i], radii[i]) < R) {
				bool found = false;
				QuadNode<index> responsibleNode = * getRoot(quad).getAppropriateLeaf(angles[i], radii[i]);

				/**
				TRACE("Getting lower bound for responsible node");
				double bound = responsibleNode.distanceLowerBound(angles[comparison], radii[comparison]);
				double actualDistance = HyperbolicSpace::getHyperbolicDistance(angles[comparison], radii[comparison], angles[i], radii[i]);
				EXPECT_GE(actualDistance, bound);
				EXPECT_TRUE(responsibleNode.responsible(angles[i], radii[i]));
	*/
				for (index j = 0; j < closeToOne.size(); j++) {
					if (closeToOne[j] == i) {
						found = true;
						break;
					}
				}
				EXPECT_TRUE(found) << "dist(" << i << "," << comparison << ") = "
						<< HyperbolicSpace::getHyperbolicDistance(angles[comparison], radii[comparison], angles[i], radii[i]) << " < " << R;
				if (!found) {
					notfound++;
					DEBUG("angle: ", angles[i], ", radius: ", radii[i], ", leftAngle: ", responsibleNode.getLeftAngle(),
							", rightAngle: ", responsibleNode.getRightAngle(), ", minR: ", responsibleNode.getMinR(), ", maxR:", responsibleNode.getMaxR());
					//DEBUG("euclidean Distance from circle center: ", center.distance(HyperbolicSpace::polarToCartesian(angles[i], radii[i])));
					DEBUG("dist(", comparison, ", leftMin)=", HyperbolicSpace::getHyperbolicDistance(angles[comparison], radii[comparison], responsibleNode.getLeftAngle(), responsibleNode.getMinR()));
					DEBUG("dist(", comparison, ", leftMax)=", HyperbolicSpace::getHyperbolicDistance(angles[comparison], radii[comparison], responsibleNode.getLeftAngle(), responsibleNode.getMaxR()));
					DEBUG("dist(", comparison, ", rightMin)=", HyperbolicSpace::getHyperbolicDistance(angles[comparison], radii[comparison], responsibleNode.getRightAngle(), responsibleNode.getMinR()));
					DEBUG("dist(", comparison, ", rightMax)=", HyperbolicSpace::getHyperbolicDistance(angles[comparison], radii[comparison], responsibleNode.getRightAngle(), responsibleNode.getMaxR()));
				}
				else {
					didfind++;
				}
			}
		}
		if (notfound > 0) {
			DEBUG("Found only ", didfind, " of ", didfind + notfound, " neighbours");
		}
	}
}

TEST_F(QuadTreeTest, testQuadTreeQuery) {
	count n = 10;
	double R = acosh((double)n/(2*M_PI)+1);
	vector<double> angles(n);
	vector<double> radii(n);
	vector<double> indices(n);
	HyperbolicSpace::fillPoints(&angles, &radii, 1, 1);
	double max = 0;
	for (index i = 0; i < n; i++) {
		indices[i] = i;
		if (radii[i] > max) {
			max = radii[i];
		}
	}

	Quadtree<index> quad(HyperbolicSpace::hyperbolicRadiusToEuclidean(R));

	for (index i = 0; i < n; i++) {
		EXPECT_GE(angles[i], 0);
		EXPECT_LT(angles[i], 2*M_PI);
		EXPECT_GE(radii[i], 0);
		EXPECT_LT(radii[i], R);
		TRACE("Added (", angles[i], ",", radii[i], ")");
		quad.addContent(i, angles[i], radii[i]);
	}

	for (index testindex = 0; testindex < 100; testindex++) {
		index query = Aux::Random::integer(n);
		if (query == n) query--;
		vector<index> lastNeighbours;
		for (double threshold = 0; threshold < R; threshold += 0.01) {
			vector<index> neighbours = quad.getCloseElements(HyperbolicSpace::polarToCartesian(angles[query], radii[query]), threshold);
			EXPECT_GE(neighbours.size(), lastNeighbours.size());
			if (neighbours.size() < lastNeighbours.size()) {
				DEBUG("Previous Neighbours: ");
				for (index neighbour : lastNeighbours) {
					DEBUG(neighbour);
				}
				DEBUG("Current Neighbours: ");
				for (index neighbour : neighbours) {
					DEBUG(neighbour);
				}
			}
			lastNeighbours = neighbours;
		}
	}
}

TEST_F(QuadTreeTest, testQuadTreeDeletion) {
	count n = 1000;
	double R = acosh((double)n/(2*M_PI)+1);
	vector<double> angles(n);
	vector<double> radii(n);
	vector<double> indices(n);
	HyperbolicSpace::fillPoints(&angles, &radii, 1, 1);
	double max = 0;
	for (index i = 0; i < n; i++) {
		indices[i] = i;
		if (radii[i] > max) {
			max = radii[i];
		}
	}

	Quadtree<index> quad(HyperbolicSpace::hyperbolicRadiusToEuclidean(R));

	for (index i = 0; i < n; i++) {
		EXPECT_GE(angles[i], 0);
		EXPECT_LT(angles[i], 2*M_PI);
		EXPECT_GE(radii[i], 0);
		EXPECT_LT(radii[i], R);
		TRACE("Added (", angles[i], ",", radii[i], ")");
		quad.addContent(i, angles[i], radii[i]);
	}

	while(indices.size() > 0) {
		index toRemove = Aux::Random::integer(indices.size());
		if (toRemove == indices.size()) toRemove--;
		assert(toRemove < indices.size());
		assert(toRemove >= 0);
		EXPECT_EQ(quad.size(), indices.size());
		bool removed = quad.removeContent(indices[toRemove], angles[toRemove], radii[toRemove]);
		EXPECT_TRUE(removed);
		EXPECT_EQ(quad.size(), indices.size()-1);
		bool removedTwice = quad.removeContent(indices[toRemove], angles[toRemove], radii[toRemove]);
		EXPECT_FALSE(removedTwice);
		EXPECT_EQ(quad.size(), indices.size()-1);
		indices.erase(indices.begin()+toRemove);
		angles.erase(angles.begin()+toRemove);
		radii.erase(radii.begin()+toRemove);
	}

	QuadNode<index> root = getRoot(quad);
	EXPECT_EQ(getChildren(root).size(), 0);//root is leaf node, coarsening worked.
}

TEST_F(QuadTreeTest, testEuclideanCircle) {
		count n = 1000;
		double R = 1;
		vector<double> angles(n);
		vector<double> radii(n);
		HyperbolicSpace::fillPoints(&angles, &radii, 1, 1);
		double max = 0;
		for (index i = 0; i < n; i++) {
			if (radii[i] > max) {
				max = radii[i];
			}
		}
		Quadtree<index> quad(max+(1-max)/4);

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
		QuadNode<index> root = getRoot(quad);
		for (index i = 0; i < 100; i++) {
			index comparison = Aux::Random::integer(n);
			Point2D<double> origin;
			Point2D<double> query = HyperbolicSpace::polarToCartesian(angles[comparison], radii[comparison]);
			double radius = Aux::Random::real(1);//this may overshoot the poincaré disc, this is intentional. I want to know what happens
			double minR = query.length() - radius;
			double maxR = query.length() + radius;
			double minPhi, maxPhi, phi_c, r_c, spread;
			if (minR < 0) {
				maxR = std::max(abs(minR), maxR);
				minR = 0;
				minPhi = 0;
				maxPhi = 2*M_PI;
			} else {
				spread = asin(radius / query.length());
				HyperbolicSpace::cartesianToPolar(query, phi_c, r_c);
				minPhi = phi_c - spread;
				maxPhi = phi_c + spread;
				/**
				 * what to do if they overlap the 2\pi line? Well, have to make two separate calls and collect
				 */
			}

			/**
			 * get Elements in Circle
			 */

			vector<index> circleDenizens;

			root.getElementsInEuclideanCircle(minPhi, maxPhi, minR, maxR, query, radius, circleDenizens);
			if (minPhi < 0) {
				root.getElementsInEuclideanCircle(2*M_PI+minPhi, 2*M_PI, minR, maxR, query, radius, circleDenizens);
			}
			if (maxPhi > 2*M_PI) {
				root.getElementsInEuclideanCircle(0, maxPhi - 2*M_PI, minR, maxR, query, radius, circleDenizens);
			}

			for (index j = 0; j < n; j++) {
				Point2D<double> comp = HyperbolicSpace::polarToCartesian(angles[j], radii[j]);
				double dist = comp.distance(query);
				if (dist < radius) {
					bool found = false;
					for (index k = 0; k < circleDenizens.size(); k++) {
						if (circleDenizens[k] == j) {
							found = true;
						}
					}
					EXPECT_TRUE(found)<< "dist(" << j << "," << comparison << ") = "
							<< dist << " < " << radius;
					if (!found) {
						DEBUG("angle: ", angles[j], ", radius: ", radii[j]);
					}
				}
			}
		}
}

TEST_F(QuadTreeTest, testQuadTreeBalance) {
	count n = 100000;
	double R = acosh((double)n/(2*M_PI)+1);
	vector<double> angles(n);
	vector<double> radii(n);
	HyperbolicSpace::fillPoints(&angles, &radii, 1, 1);
	double max = 0;
	for (index i = 0; i < n; i++) {
		if (radii[i] > max) {
			max = radii[i];
		}
	}

	Quadtree<index> quad(HyperbolicSpace::hyperbolicRadiusToEuclidean(R));

	for (index i = 0; i < n; i++) {
		EXPECT_GE(angles[i], 0);
		EXPECT_LT(angles[i], 2*M_PI);
		EXPECT_GE(radii[i], 0);
		EXPECT_LT(radii[i], R);
		TRACE("Added (", angles[i], ",", radii[i], ")");
		quad.addContent(i, angles[i], radii[i]);
	}

	QuadNode<index> root = getRoot(quad);
	vector<QuadNode<index> > children = getChildren(root);
	EXPECT_EQ(children.size(), 4);

	EXPECT_LE(children[0].size(), 2*children[1].size());
	EXPECT_LE(children[0].size(), 2*children[3].size());

	EXPECT_LE(children[1].size(), 2*children[0].size());
	EXPECT_LE(children[1].size(), 2*children[2].size());

	EXPECT_LE(children[2].size(), 2*children[1].size());
	EXPECT_LE(children[2].size(), 2*children[3].size());

	EXPECT_LE(children[3].size(), 2*children[2].size());
	EXPECT_LE(children[3].size(), 2*children[0].size());

}

} /* namespace NetworKit */

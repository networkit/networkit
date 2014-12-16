/*
 * QuadTreeTest.cpp
 *
 *  Created on: 28.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#include <stack>
#include <cmath>
#include <algorithm>

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

/**
 * Test whether the elements returned by a quadtree range query are indeed those whose hyperbolic distance to the query point is below a threshold
 */
TEST_F(QuadTreeTest, testQuadTreeHyperbolicCircle) {
	count n = 1000;
	double R = 1;
	vector<double> angles(n);
	vector<double> radii(n);
	HyperbolicSpace::fillPoints(angles, radii, 1, 1);
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
	EXPECT_EQ(quad.size(), n);
	for (index testindex = 0; testindex < 100; testindex++) {
		index comparison = Aux::Random::integer(n-1);
		Point2D<double> origin;
		Point2D<double> query = HyperbolicSpace::polarToCartesian(angles[comparison], radii[comparison]);
		DEBUG("Using ", comparison, " at (", angles[comparison], ",", radii[comparison], ") as query point");

		vector<index> closeToOne = quad.getElementsInHyperbolicCircle(HyperbolicSpace::polarToCartesian(angles[comparison], radii[comparison]), R);
		EXPECT_LE(closeToOne.size(), n);

		for (index i = 0; i < closeToOne.size(); i++) {
			//no corrupt indices
			ASSERT_LE(closeToOne[i], n);

			//close results should actually be close
			EXPECT_LE(HyperbolicSpace::poincareMetric(angles[comparison], radii[comparison], angles[closeToOne[i]], radii[closeToOne[i]]), R);
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
			if (HyperbolicSpace::poincareMetric(angles[comparison], radii[comparison], angles[i], radii[i]) < R) {
				bool found = false;
				QuadNode<index> responsibleNode = getRoot(quad).getAppropriateLeaf(angles[i], radii[i]);

				for (index j = 0; j < closeToOne.size(); j++) {
					if (closeToOne[j] == i) {
						found = true;
						break;
					}
				}
				EXPECT_TRUE(found) << "dist(" << i << "," << comparison << ") = "
						<< HyperbolicSpace::poincareMetric(angles[comparison], radii[comparison], angles[i], radii[i]) << " < " << R;
				if (!found) {
					notfound++;
					DEBUG("angle: ", angles[i], ", radius: ", radii[i], ", leftAngle: ", responsibleNode.getLeftAngle(),
							", rightAngle: ", responsibleNode.getRightAngle(), ", minR: ", responsibleNode.getMinR(), ", maxR:", responsibleNode.getMaxR());
					//DEBUG("euclidean Distance from circle center: ", center.distance(HyperbolicSpace::polarToCartesian(angles[i], radii[i])));
					DEBUG("dist(", comparison, ", leftMin)=", HyperbolicSpace::poincareMetric(angles[comparison], radii[comparison], responsibleNode.getLeftAngle(), responsibleNode.getMinR()));
					DEBUG("dist(", comparison, ", leftMax)=", HyperbolicSpace::poincareMetric(angles[comparison], radii[comparison], responsibleNode.getLeftAngle(), responsibleNode.getMaxR()));
					DEBUG("dist(", comparison, ", rightMin)=", HyperbolicSpace::poincareMetric(angles[comparison], radii[comparison], responsibleNode.getRightAngle(), responsibleNode.getMinR()));
					DEBUG("dist(", comparison, ", rightMax)=", HyperbolicSpace::poincareMetric(angles[comparison], radii[comparison], responsibleNode.getRightAngle(), responsibleNode.getMaxR()));
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

/**
 * Gradually increase the distance threshold and check whether the number of neighbours increases monotonically. Necessary foundation for the dynamic hyperbolic generator.
 */
TEST_F(QuadTreeTest, testQuadTreeThresholdGrowth) {
	count n = 100;
	double R = HyperbolicSpace::hyperbolicAreaToRadius(n);
	vector<double> angles(n);
	vector<double> radii(n);
	vector<double> indices(n);
	HyperbolicSpace::fillPoints(angles, radii, 1, 1);
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
			vector<index> neighbours = quad.getElementsInHyperbolicCircle(HyperbolicSpace::polarToCartesian(angles[query], radii[query]), threshold);
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

/**
 * Insert nodes into Quadtree and successively delete all of them, check if resulting tree is empty
 */
TEST_F(QuadTreeTest, testQuadTreeDeletion) {
	count n = 1000;
	double R = HyperbolicSpace::hyperbolicAreaToRadius(n);
	vector<double> angles(n);
	vector<double> radii(n);
	vector<double> indices(n);
	HyperbolicSpace::fillPoints(angles, radii, 1, 1);
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
		//pick random point which is not yet deleted
		index toRemove = Aux::Random::integer(indices.size());
		if (toRemove == indices.size()) toRemove--;
		assert(toRemove < indices.size());
		assert(toRemove >= 0);

		//remove content at point
		EXPECT_EQ(quad.size(), indices.size());
		bool removed = quad.removeContent(indices[toRemove], angles[toRemove], radii[toRemove]);
		EXPECT_TRUE(removed);
		EXPECT_EQ(quad.size(), indices.size()-1);

		//removing twice should not work
		bool removedTwice = quad.removeContent(indices[toRemove], angles[toRemove], radii[toRemove]);
		EXPECT_FALSE(removedTwice);
		EXPECT_EQ(quad.size(), indices.size()-1);

		//mark point as removed in query list
		indices.erase(indices.begin()+toRemove);
		angles.erase(angles.begin()+toRemove);
		radii.erase(radii.begin()+toRemove);
	}

	QuadNode<index> root = getRoot(quad);
	//if root is leaf node, the coarsening worked
	EXPECT_EQ(getChildren(root).size(), 0);
}

/**
 * Test whether the points found by a Euclidean range query on the quadtree root are exactly those whose Euclidean distance to the query point is smaller than the threshold.
 */
TEST_F(QuadTreeTest, testEuclideanCircle) {
	count n = 1000;
	double R = 1;
	vector<double> angles(n);
	vector<double> radii(n);
	HyperbolicSpace::fillPoints(angles, radii, 1, 1);
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

		root.getElementsInEuclideanCircle(query, radius, circleDenizens, minPhi, maxPhi, minR, maxR);
		if (minPhi < 0) {
			root.getElementsInEuclideanCircle(query, radius, circleDenizens, 2*M_PI+minPhi, 2*M_PI, minR, maxR);
		}
		if (maxPhi > 2*M_PI) {
			root.getElementsInEuclideanCircle(query, radius, circleDenizens, 0, maxPhi - 2*M_PI, minR, maxR);
		}

		//check whether bounds were correct by calling again without bounds and comparing
		vector<index> alternateDenizens;
		root.getElementsInEuclideanCircle(query, radius, alternateDenizens);

		EXPECT_EQ(circleDenizens.size(), alternateDenizens.size());


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

/**
 * Test whether the theoretical splitting rules for different point distributions succeed in creating a balanced tree
 */
TEST_F(QuadTreeTest, testQuadTreeBalance) {
	count n = 100000;
	double s =1;
	double alpha = 1;
	double R = s*HyperbolicSpace::hyperbolicAreaToRadius(n);
	vector<double> angles(n);
	vector<double> radii(n);
	HyperbolicSpace::fillPoints(angles, radii, s, alpha);
	double max = 0;
	for (index i = 0; i < n; i++) {
		if (radii[i] > max) {
			max = radii[i];
		}
	}

	Quadtree<index> quad(HyperbolicSpace::hyperbolicRadiusToEuclidean(R),true,alpha);

	for (index i = 0; i < n; i++) {
		EXPECT_GE(angles[i], 0);
		EXPECT_LT(angles[i], 2*M_PI);
		EXPECT_GE(radii[i], 0);
		EXPECT_LT(radii[i], R);
		TRACE("Added (", angles[i], ",", radii[i], ")");
		quad.addContent(i, angles[i], radii[i]);
	}

	QuadNode<index> root = getRoot(quad);

	//visit tree
	std::stack<QuadNode<index> > toAnalyze;
	toAnalyze.push(root);
	while (!toAnalyze.empty()) {
		QuadNode<index> current = toAnalyze.top();
		toAnalyze.pop();
		if (current.height() > 1) {
			vector<QuadNode<index> > children = getChildren(current);
			EXPECT_EQ(children.size(), 4);

			EXPECT_LE(children[0].size(), 2*children[1].size());
			EXPECT_LE(children[0].size(), 2*children[3].size());

			EXPECT_LE(children[1].size(), 2*children[0].size());
			EXPECT_LE(children[1].size(), 2*children[2].size());

			EXPECT_LE(children[2].size(), 2*children[1].size());
			EXPECT_LE(children[2].size(), 2*children[3].size());

			EXPECT_LE(children[3].size(), 2*children[2].size());
			EXPECT_LE(children[3].size(), 2*children[0].size());

			EXPECT_LE(children[0].height(), children[1].height()+1);
			EXPECT_LE(children[0].height(), children[3].height()+1);

			EXPECT_LE(children[1].height(), children[0].height()+1);
			EXPECT_LE(children[1].height(), children[2].height()+1);

			EXPECT_LE(children[2].height(), children[1].height()+1);
			EXPECT_LE(children[2].height(), children[3].height()+1);

			EXPECT_LE(children[3].height(), children[2].height()+1);
			EXPECT_LE(children[3].height(), children[0].height()+1);
			for (auto child : children) {
				toAnalyze.push(child);
			}
		}
	}
}

/**
 * No testing yet, just debug output
 */
TEST_F(QuadTreeTest, tryQuadTreeCutLeaves) {
	count n = 1000000;
	count trials = 20;
	double s =1;
	double alpha = 1;
	double t = 1;
	double R = s*HyperbolicSpace::hyperbolicAreaToRadius(n);
	double threshold = t*R;
	vector<double> angles(n);
	vector<double> radii(n);
	HyperbolicSpace::fillPoints(angles, radii, s, alpha);

	for (index capexp = 1; capexp < log(n)/log(4); capexp++) {
		count capacity = pow(4,capexp);
		Quadtree<index> quad(HyperbolicSpace::hyperbolicRadiusToEuclidean(R),false,alpha,capacity,true);

		for (index i = 0; i < n; i++) {
			quad.addContent(i, angles[i], radii[i]);
		}

		count sumIncluded = 0;
		count sumCut = 0;
		count sumUNcomp = 0;
		count sumNcomp = 0;
		count totalEdges = 0;

		for (index e = 0; e < trials; e++) {
			index q = Aux::Random::integer(n);
			quad.resetCounter();
			vector<index> neighbours = quad.getElementsInHyperbolicCircle(HyperbolicSpace::polarToCartesian(angles[q], radii[q]), threshold);
			QuadNode<index> root = getRoot(quad);
			count included = root.countIncluded();
			count cut = root.countCut();
			count uncomp = root.countUnnecessaryComparisonsInCutLeaves();
			count ncomp = root.countNecessaryComparisonsInCutLeaves();
			/*
			DEBUG("Node: ", q);
			DEBUG("Degree: ", neighbours.size());
			DEBUG("Included leaves: ", included);
			DEBUG("Cut leaves: ", cut);
			DEBUG("Unnecessary comparisons in cut leaves: ", uncomp);
			DEBUG("Necessary comparisons in cut leaves: ", ncomp);
			*/
			sumIncluded += included;
			sumCut += cut;
			sumUNcomp += uncomp;
			sumNcomp += ncomp;
			totalEdges += neighbours.size();
		}
		DEBUG("Capacity:", capacity);
		DEBUG("Number of Leaves:", quad.countLeaves(), "(", pow(4,ceil(log(n/capacity)/log(4))), ")");
		DEBUG("Average included leaves: ", sumIncluded/trials);
		DEBUG("Average cut leaves: ", sumCut/trials);
		DEBUG("Avg unnecessary comparisons in cut leaves: ", sumUNcomp/trials);
		DEBUG("Avg necessary comparisons in cut leaves: ", sumNcomp/trials);
		DEBUG("Avg Total edges: ", totalEdges/trials);
	}
}

/**
 * No testing yet, just debug output
 */
TEST_F(QuadTreeTest, testQuadTreeCutLeaves) {
	count n = 100000;
	count capacity = 1000;
	count trials = 20;
	double s =1;
	double alpha = 1;
	double t = 1;
	double R = s*HyperbolicSpace::hyperbolicAreaToRadius(n);
	double threshold = t*R;
	vector<double> angles(n);
	vector<double> radii(n);
	HyperbolicSpace::fillPoints(angles, radii, s, alpha);

	Quadtree<index> quad(HyperbolicSpace::hyperbolicRadiusToEuclidean(R),false,alpha,capacity,true);

	for (index i = 0; i < n; i++) {
		quad.addContent(i, angles[i], radii[i]);
	}

	count sumIncluded = 0;
	count sumCut = 0;
	count sumUNcomp = 0;
	count sumNcomp = 0;
	count totalEdges = 0;

	for (index e = 0; e < trials; e++) {
		index q = Aux::Random::integer(n);
		quad.resetCounter();
		vector<index> neighbours = quad.getElementsInHyperbolicCircle(HyperbolicSpace::polarToCartesian(angles[q], radii[q]), threshold);
		QuadNode<index> root = getRoot(quad);
		count included = root.countIncluded();
		count cut = root.countCut();
		count uncomp = root.countUnnecessaryComparisonsInCutLeaves();
		count ncomp = root.countNecessaryComparisonsInCutLeaves();
		DEBUG("Node: ", q);
		DEBUG("Degree: ", neighbours.size());
		DEBUG("Included leaves: ", included);
		DEBUG("Cut leaves: ", cut);
		DEBUG("Unnecessary comparisons in cut leaves: ", uncomp);
		DEBUG("Necessary comparisons in cut leaves: ", ncomp);
		sumIncluded += included;
		sumCut += cut;
		sumUNcomp += uncomp;
		sumNcomp += ncomp;
		totalEdges += neighbours.size();
	}
	DEBUG("Number of Leaves:", quad.countLeaves());
	DEBUG("Total included leaves: ", sumIncluded);
	DEBUG("Total cut leaves: ", sumCut);
	DEBUG("Unnecessary comparisons in cut leaves: ", sumUNcomp);
	DEBUG("Necessary comparisons in cut leaves: ", sumNcomp);
	DEBUG("Total edges: ", totalEdges);
}

TEST_F(QuadTreeTest, testParallelQuadtreeConstruction) {
	count n = 1000000;
	double s = 1;
	Quadtree<index> quad(n,s);
	EXPECT_EQ(quad.size(), n);
	EXPECT_GE(quad.height(), log(n/1000)/log(4));
	EXPECT_GE(quad.countLeaves(), n/1000);
	quad.reindex();
	vector<index> elements = quad.getElements();
	EXPECT_EQ(elements.size(), n);
	for (index i = 0; i < elements.size(); i++) {
		EXPECT_EQ(elements[i], i);
	}
}

TEST_F(QuadTreeTest, testSequentialQuadtreeConstruction) {
	count n = 1000000;
	count capacity = 1000;
	double s =1;
	double alpha = 1;
	double R = s*HyperbolicSpace::hyperbolicAreaToRadius(n);
	vector<double> angles(n);
	vector<double> radii(n);
	HyperbolicSpace::fillPoints(angles, radii, s, alpha);

	Quadtree<index> quad(HyperbolicSpace::hyperbolicRadiusToEuclidean(R),false,alpha,capacity,true);

	for (index i = 0; i < n; i++) {
		quad.addContent(i, angles[i], radii[i]);
	}
	EXPECT_EQ(quad.size(), n);

	quad.trim();
	quad.sortPointsInLeaves();
	vector<double> anglecopy;
	vector<double> radiicopy;
	quad.extractCoordinates(anglecopy, radiicopy);
	EXPECT_EQ(anglecopy.size(), n);
	EXPECT_EQ(radiicopy.size(), n);
	//EXPECT_TRUE(std::is_permutation(angles.begin(), angles.end(), anglecopy.begin()));//too slow!
	//EXPECT_TRUE(std::is_permutation(radii.begin(), radii.end(), radiicopy.begin()));
}


} /* namespace NetworKit */

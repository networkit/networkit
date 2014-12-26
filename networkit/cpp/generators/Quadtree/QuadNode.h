/*
 * QuadNode.h
 *
 *  Created on: 21.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#ifndef QUADNODE_H_
#define QUADNODE_H_

#include <vector>
#include <algorithm>
#include <assert.h>
#include "../../auxiliary/Log.h"
#include "../../geometric/HyperbolicSpace.h"

using std::vector;
using std::min;
using std::cos;

namespace NetworKit {

template <class T>
class QuadNode {
	friend class QuadTreeTest;
private:
	double leftAngle;
	double minR;
	double rightAngle;
	double maxR;
	Point2D<double> a,b,c,d;
	unsigned capacity;
	static const unsigned coarsenLimit = 4;
	double minRegion;//the minimal region a QuadNode should cover. If it is smaller, don't bother splitting up.
	count elements;
	std::vector<T> content;
	std::vector<Point2D<double> > positions;
	std::vector<double> angles;
	std::vector<double> radii;
	bool isLeaf;
	bool splitTheoretical;
	double alpha;
	double balance;
	index ID;
	double lowerBoundR;

public:
	std::vector<QuadNode> children;

	QuadNode() {
		//This should never be called.
		leftAngle = 0;
		rightAngle = 0;
		minR = 0;
		maxR = 0;
		capacity = 20;
		isLeaf = true;
		minRegion = 0;
		elements = 0;
		balance = 0.5;
		alpha = 1;
		lowerBoundR = maxR;
	}

	/**
	 * Construct a QuadNode for polar coordinates.
	 *
	 *
	 * @param leftAngle Minimal angular coordinate of region, in radians from 0 to 2\pi
	 * @param rightAngle Maximal angular coordinate of region, in radians from 0 to 2\pi
	 * @param minR Minimal radial coordinate of region, between 0 and 1
	 * @param maxR Maximal radial coordinate of region, between 0 and 1
	 * @param capacity Number of points a leaf cell can store before splitting
	 * @param minDiameter Minimal diameter of a quadtree node. If the node is already smaller, don't split even if over capacity. Default is 0
	 * @param splitTheoretical Whether to split in a theoretically optimal way or in a way to decrease measured running times
	 * @param alpha dispersion Parameter of the point distribution. Only has an effect if theoretical split is true
	 * @param diagnostics Count how many necessary and unnecessary comparisons happen in leaf cells? Will cause race condition and false sharing in parallel use
	 *
	 */
	QuadNode(double leftAngle, double minR, double rightAngle, double maxR, unsigned capacity, double minDiameter, bool splitTheoretical = false, double alpha = 1, double balance = 0.5) {
		this->leftAngle = leftAngle;
		this->minR = minR;
		this->maxR = maxR;
		this->rightAngle = rightAngle;
		this->a = HyperbolicSpace::polarToCartesian(leftAngle, minR);
		this->b = HyperbolicSpace::polarToCartesian(rightAngle, minR);
		this->c = HyperbolicSpace::polarToCartesian(rightAngle, maxR);
		this->d = HyperbolicSpace::polarToCartesian(leftAngle, maxR);
		this->capacity = capacity;
		this->minRegion = minDiameter;
		this->alpha = alpha;
		this->splitTheoretical = splitTheoretical;
		this->balance = balance;
		assert(balance > 0);
		assert(balance < 1);
		this->lowerBoundR = maxR;
		isLeaf = true;
		elements = 0;
	}

	void split() {
		assert(isLeaf);
		//heavy lifting: split up!
		double middleAngle = (rightAngle - leftAngle) / 2 + leftAngle;
		/**
		 * we want to make sure the space is evenly divided to obtain a balanced tree
		 * Simply halving the radius will cause a larger space for the outer Quadnode, resulting in an unbalanced tree
		 */

		double middleR;
		if (splitTheoretical) {
			double hyperbolicOuter = HyperbolicSpace::EuclideanRadiusToHyperbolic(maxR);
			double hyperbolicInner = HyperbolicSpace::EuclideanRadiusToHyperbolic(minR);
			double hyperbolicMiddle = acosh((1-balance)*cosh(alpha*hyperbolicOuter) + balance*cosh(alpha*hyperbolicInner))/alpha;
			middleR = HyperbolicSpace::hyperbolicRadiusToEuclidean(hyperbolicMiddle);
		} else {
			double nom = maxR - minR;
			double denom = pow((1-maxR*maxR)/(1-minR*minR), 0.5)+1;
			middleR = nom/denom + minR;
		}

		//one could also use the median here. Results in worse asymptotical complexity, but maybe better runtime?

		assert(middleR < maxR);
		assert(middleR > minR);

		QuadNode southwest(leftAngle, minR, middleAngle, middleR, capacity, minRegion, splitTheoretical, alpha, balance);
		QuadNode southeast(middleAngle, minR, rightAngle, middleR, capacity, minRegion, splitTheoretical, alpha, balance);
		QuadNode northwest(leftAngle, middleR, middleAngle, maxR, capacity, minRegion, splitTheoretical, alpha, balance);
		QuadNode northeast(middleAngle, middleR, rightAngle, maxR, capacity, minRegion, splitTheoretical, alpha, balance);
		children = {southwest, southeast, northwest, northeast};
		isLeaf = false;
	}

	/**
	 * Add a point at polar coordinates (angle, R) with content input. May split node if capacity is full
	 *
	 * @param input arbitrary content, in our case an index
	 * @param angle angular coordinate of point, between 0 and 2 pi.
	 * @param R radial coordinate of point, between 0 and 1.
	 */
	void addContent(T input, double angle, double R) {
		assert(this->responsible(angle, R));
		if (lowerBoundR > R) lowerBoundR = R;
		if (isLeaf) {
			if (content.size() + 1 < capacity ||  HyperbolicSpace::poincareMetric(leftAngle, minR, rightAngle, maxR) < minRegion) {
				content.push_back(input);
				angles.push_back(angle);
				radii.push_back(R);
				Point2D<double> pos = HyperbolicSpace::polarToCartesian(angle, R);
				positions.push_back(pos);
			} else {

				split();

				for (uint i = 0; i < content.size(); i++) {
					this->addContent(content[i], angles[i], radii[i]);
				}

				content.clear();
				this->addContent(input, angle, R);
			}
		}
		else {
			assert(children.size() > 0);
			for (uint i = 0; i < children.size(); i++) {
				if (children[i].responsible(angle, R)) {
					children[i].addContent(input, angle, R);
					break;
				}
			}
		}
		elements++;
	}

	/**
	 * Remove content at polar coordinates (angle, R). May cause coarsening of the quadtree
	 *
	 * @param input Content to be removed
	 * @param angle Angular coordinate
	 * @param R Radial coordinate
	 *
	 * @return True if content was found and removed, false otherwise
	 */
	bool removeContent(T input, double angle, double R) {
		if (!responsible(angle, R)) return false;
		if (isLeaf) {
			index i = 0;
			for (; i < content.size(); i++) {
				if (content[i] == input) break;
			}
			if (i < content.size()) {
				assert(angles[i] == angle);
				assert(radii[i] == R);
				//remove element
				content.erase(content.begin()+i);
				positions.erase(positions.begin()+i);
				angles.erase(angles.begin()+i);
				radii.erase(radii.begin()+i);
				return true;
			} else {
				return false;
			}
		}
		else {
			bool removed = false;
			bool allLeaves = true;
			assert(children.size() > 0);
			for (index i = 0; i < children.size(); i++) {
				if (!children[i].isLeaf) allLeaves = false;
				if (children[i].removeContent(input, angle, R)) {
					assert(!removed);
					removed = true;
				}
			}
			//coarsen?
			if (removed && allLeaves && size() < coarsenLimit) {
				//coarsen!!
				//why not assert empty containers and then insert directly?
				vector<T> allContent;
				vector<Point2D<double> > allPositions;
				vector<double> allAngles;
				vector<double> allRadii;
				for (index i = 0; i < children.size(); i++) {
					allContent.insert(allContent.end(), children[i].content.begin(), children[i].content.end());
					allPositions.insert(allPositions.end(), children[i].positions.begin(), children[i].positions.end());
					allAngles.insert(allAngles.end(), children[i].angles.begin(), children[i].angles.end());
					allRadii.insert(allRadii.end(), children[i].radii.begin(), children[i].radii.end());
				}
				assert(allContent.size() == allPositions.size());
				assert(allContent.size() == allAngles.size());
				assert(allContent.size() == allRadii.size());
				children.clear();
				content.swap(allContent);
				positions.swap(allPositions);
				angles.swap(allAngles);
				radii.swap(allRadii);
				isLeaf = true;
			}

			return removed;
		}
	}


	/**
	 * Check whether the region managed by this node lies outside of an Euclidean circle.
	 *
	 * @param query Center of the Euclidean query circle, given in Cartesian coordinates
	 * @param radius Radius of the Euclidean query circle
	 *
	 * @return True if the region managed by this node lies completely outside of the circle
	 */
	bool outOfReach(Point2D<double> query, double radius) const {
		double phi, r;
		HyperbolicSpace::cartesianToPolar(query, phi, r);
		if (responsible(phi, r)) return false;

		//get four edge points
		double topDistance, bottomDistance, leftDistance, rightDistance;

		if (phi < leftAngle || phi > rightAngle) {
			topDistance = min(c.distance(query), d.distance(query));
		} else {
			topDistance = abs(r - maxR);
		}
		if (topDistance <= radius) return false;
		if (phi < leftAngle || phi > rightAngle) {
			bottomDistance = min(a.distance(query), b.distance(query));
		} else {
			bottomDistance = abs(r - minR);
		}
		if (bottomDistance <= radius) return false;

		double minDistanceR = r*cos(abs(phi-leftAngle));
		if (minDistanceR > minR && minDistanceR < maxR) {
			leftDistance = query.distance(HyperbolicSpace::polarToCartesian(phi, minDistanceR));
		} else {
			leftDistance = min(a.distance(query), d.distance(query));
		}
		if (leftDistance <= radius) return false;

		minDistanceR = r*cos(abs(phi-rightAngle));
		if (minDistanceR > minR && minDistanceR < maxR) {
			rightDistance = query.distance(HyperbolicSpace::polarToCartesian(phi, minDistanceR));
		} else {
			rightDistance = min(b.distance(query), c.distance(query));
		}
		if (rightDistance <= radius) return false;
		return true;
	}

	/**
	 * Check whether the region managed by this node lies outside of an Euclidean circle.
	 * Functionality is the same as in the method above, but it takes polar coordinates instead of Cartesian ones
	 *
	 * @param angle_c Angular coordinate of the Euclidean query circle's center
	 * @param r_c Radial coordinate of the Euclidean query circle's center
	 * @param radius Radius of the Euclidean query circle
	 *
	 * @return True if the region managed by this node lies completely outside of the circle
	 */
	bool outOfReach(double angle_c, double r_c, double radius) const {
		if (responsible(angle_c, r_c)) return false;
		Point2D<double> query = HyperbolicSpace::polarToCartesian(angle_c, r_c);
		return outOfReach(query, radius);
	}

	/**
	 * Does the point at (angle, r) fall inside the region managed by this QuadNode?
	 *
	 * @param angle Angular coordinate of input point
	 * @param r Radial coordinate of input points
	 *
	 * @return True if input point lies within the region of this QuadNode
	 */
	bool responsible(double angle, double r) const {
		return (angle >= leftAngle && angle < rightAngle && r >= minR && r < maxR);
	}

	/**
	 * Get all Elements in this QuadNode or a descendant of it
	 *
	 * @return vector of content type T
	 */
	std::vector<T> getElements() const {
		if (isLeaf) {
			return content;
		} else {
			vector<T> result;
			for (uint i = 0; i < children.size(); i++) {
				std::vector<T> subresult = children[i].getElements();
				result.insert(result.end(), subresult.begin(), subresult.end());
			}
			return result;
		}
	}

	void getCoordinates(vector<double> &anglesContainer, vector<double> &radiiContainer) const {
		assert(angles.size() == radii.size());
		if (isLeaf) {
			anglesContainer.insert(anglesContainer.end(), angles.begin(), angles.end());
			radiiContainer.insert(radiiContainer.end(), radii.begin(), radii.end());
		}
		else {
			for (uint i = 0; i < children.size(); i++) {
				children[i].getCoordinates(anglesContainer, radiiContainer);
			}
		}
	}

	/**
	 * Don't use this!
	 * Code is still in here for a unit test.
	 *
	 * Get copy of the leaf cell responsible for a point at (angle, r).
	 * Expensive because it copies the whole subtree, causes assertion failure if called with the wrong arguments
	 *
	 * @param angle Angular coordinate of point
	 * @param r Radial coordinate of point
	 *
	 * @return Copy of leaf cell containing point, or dummy cell not responsible for point
	 *
	 */
	QuadNode<T> getAppropriateLeaf(double angle, double r) {
		assert(this->responsible(angle, r));
		if (isLeaf) return *this;
		else {
			for (uint i = 0; i < children.size(); i++) {
				bool foundResponsibleChild = false;
				if (children[i].responsible(angle, r)) {
					assert(foundResponsibleChild == false);
					foundResponsibleChild = true;
					return children[i].getAppropriateLeaf(angle, r);
				}
			}
			DEBUG("No responsible child for (", angle, ", ", r, ") found.");
			assert(false);
			//to make compiler happy:
			QuadNode<T> dummy;
			return dummy;
		}
	}

	/**
	 * Main query method, get points lying in a Euclidean circle around the center point.
	 * Optional limits can be given to get a different result or to reduce unnecessary comparisons
	 *
	 * Elements are pushed onto a vector which is a required argument. This is done to reduce copying
	 *
	 * Safe to call in parallel if diagnostics are disabled
	 *
	 * @param center Center of the query circle
	 * @param radius Radius of the query circle
	 * @param result Reference to the vector where the results will be stored
	 * @param minAngle Optional value for the minimum angular coordinate of the query region
	 * @param maxAngle Optional value for the maximum angular coordinate of the query region
	 * @param lowR Optional value for the minimum radial coordinate of the query region
	 * @param highR Optional value for the maximum radial coordinate of the query region
	 */
	void getElementsInEuclideanCircle(Point2D<double> center, double radius, vector<T> &result, double minAngle=0, double maxAngle=2*M_PI, double lowR=0, double highR = 1) {
		if (minAngle >= rightAngle || maxAngle <= leftAngle || lowR >= maxR || highR < lowerBoundR) return;
		if (outOfReach(center, radius)) {
			return;
		}
		double phi_c, r_c;
		HyperbolicSpace::cartesianToPolar(center, phi_c, r_c);
		highR = min(highR, HyperbolicSpace::maxRinSlice(leftAngle, rightAngle, phi_c, r_c, radius));
		if (highR < lowerBoundR) return;

		if (isLeaf) {
			const double rsq = radius*radius;
			const double queryX = center[0];
			const double queryY = center[1];
			const count cSize = content.size();

			for (int i = 0; i < cSize; i++) {
				const double deltaX = positions[i][0] - queryX;
				const double deltaY = positions[i][1] - queryY;
				if (deltaX*deltaX + deltaY*deltaY < rsq) {
					result.push_back(content[i]);
				}
			}
		}	else {
			for (uint i = 0; i < children.size(); i++) {
				children[i].getElementsInEuclideanCircle(center, radius, result, minAngle, maxAngle, lowR, highR);
			}
		}
	}

	/**
	 * Shrink all vectors in this subtree to fit the content.
	 * Call after quadtree construction is complete, causes better memory usage and cache efficiency
	 */
	void trim() {
		content.shrink_to_fit();
		positions.shrink_to_fit();
		angles.shrink_to_fit();
		radii.shrink_to_fit();
		if (!isLeaf) {
			for (index i = 0; i < children.size(); i++) {
				children[i].trim();
			}
		}
	}

	/**
	 * Number of points lying in the region managed by this QuadNode
	 */
	count size() const {
		count result = isLeaf ? content.size() : 0;
		for (auto child : children) result += child.size();
		return result;
	}

	/**
	 * Height of subtree hanging from this QuadNode
	 */
	count height() const {
		count result = 1;//if leaf node, the children loop will not execute
		for (auto child : children) result = std::max(result, child.height()+1);
		return result;
	}

	/**
	 * Leaf cells in the subtree hanging from this QuadNode
	 */
	count countLeaves() const {
		if (isLeaf) return 1;
		count result = 0;
		for (index i = 0; i < children.size(); i++) {
			result += children[i].countLeaves();
		}
		return result;
	}

	double getLeftAngle() const {
		return leftAngle;
	}

	double getRightAngle() const {
		return rightAngle;
	}

	double getMinR() const {
		return minR;
	}

	double getMaxR() const {
		return maxR;
	}

	void setID(index id) {
		this->ID = id;
	}

	index getID() const {
		return ID;
	}

	index indexSubtree(index nextID) {
		if (isLeaf) {
			this->ID = nextID;
			return nextID +1;
		}
		index result = nextID;
		for (int i = 0; i < 4; i++) {
			result = children[i].indexSubtree(result);
		}
		return result;
	}

	index getCellID(double phi, double r) {
		if (!responsible(phi, r)) return -1;
		if (isLeaf) return getID();
		else {
			for (int i = 0; i < 4; i++) {
				index childresult = children[i].getCellID(phi, r);
				if (childresult >= 0) return childresult;
			}
			assert(false); //if responsible
			return -1;
		}
	}

	index getMaxIDInSubtree() {
		if (isLeaf) return getID();
		else {
			index result = -1;
			for (int i = 0; i < 4; i++) {
				result = std::max(children[i].getMaxIDInSubtree(), result);
			}
			return result;
		}
	}

	count reindex(count offset) {
		if (isLeaf)
		{
			#pragma omp task
			{
				index p = offset;
				std::generate(content.begin(), content.end(), [&p](){return p++;});
			}
			offset += size();
		} else {
			for (int i = 0; i < 4; i++) {
				offset = children[i].reindex(offset);
			}
		}
		return offset;
	}

	void sortPointsInLeaves() {
		if (isLeaf) {
			#pragma omp task
			{
				count cs = content.size();
				vector<index> permutation(cs);

				index p = 0;
				std::generate(permutation.begin(), permutation.end(), [&p](){return p++;});

				//can probably be parallelized easily, but doesn't bring much benefit
				std::sort(permutation.begin(), permutation.end(), [this](index i, index j){return angles[i] < angles[j];});

				//There ought to be a way to do this more elegant with some algorithm header, but I didn't find any

				std::vector<T> contentcopy(cs);
				std::vector<Point2D<double> > positioncopy(cs);
				std::vector<double> anglecopy(cs);
				std::vector<double> radiicopy(cs);

				for (index i = 0; i < cs; i++) {
					const index perm = permutation[i];
					contentcopy[i] = content[perm];
					positioncopy[i] = positions[perm];
					anglecopy[i] = angles[perm];
					radiicopy[i] = radii[perm];
				}

				content.swap(contentcopy);
				positions.swap(positioncopy);
				angles.swap(anglecopy);
				radii.swap(radiicopy);
			}

		} else {
			for (int i = 0; i < 4; i++) {
				children[i].sortPointsInLeaves();
			}
		}
	}
};
}

#endif /* QUADNODE_H_ */

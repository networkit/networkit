/*
 * Octree.h
 *
 *  Created on: Apr 21, 2016
 *      Author: Michael Wegner
 */

#ifndef NETWORKIT_CPP_VIZ_OCTREE_H_
#define NETWORKIT_CPP_VIZ_OCTREE_H_

#include <vector>

#include "Point.h"
#include "../algebraic/Vector.h"

#include "../auxiliary/Log.h"

namespace NetworKit {

/**
 * Bounding box used by the Octree class.
 */
template<typename T>
struct BoundingBox {
public:
	BoundingBox() : center(Point<T>()), sideLength(0), halfSideLength(0), sqSideLength(0), dimension(0) {}
	BoundingBox(const Point<T>& center, const T sideLength) : center(center), sideLength(sideLength), halfSideLength(sideLength/2.0), sqSideLength(sideLength*sideLength), dimension(center.getDimensions()) {}
	BoundingBox(const BoundingBox<T>& other) : center(other.center), sideLength(other.sideLength), halfSideLength(other.halfSideLength), sqSideLength(other.sqSideLength), dimension(other.dimension) {}

	void setCenter(const Point<T>& center) {
		this->center = center;
		dimension = center.getDimensions();
	}

	inline Point<T>& getCenter() {
		return center;
	}

	void setSideLength(T sideLength) {
		this->sideLength = sideLength;
		this->halfSideLength = sideLength/2.0;
		this->sqSideLength = sideLength * sideLength;
	}

	inline T getSideLength() const {
		return sideLength;
	}

	inline T getHalfSideLength() const {
		return halfSideLength;
	}

	inline T getSqSideLength() const {
		return sqSideLength;
	}

	bool contains(const Point<T>& point) const {
		for (index d = 0; d < dimension; ++d) {
			if (center[d] - halfSideLength > point[d] || point[d] > center[d] + halfSideLength) {
				return false;
			}
		}

		return true;
	}

private:
	Point<T> center;
	T sideLength;
	T halfSideLength;
	T sqSideLength;
	count dimension;
};

/**
 * Node in the octree data structure.
 */
template<typename T>
struct OctreeNode {
	count weight;
	Point<T> centerOfMass;
	std::vector<OctreeNode> children;
	BoundingBox<T> bBox;

	OctreeNode() : weight(0), centerOfMass({0,0}), children({}), bBox(BoundingBox<T>()) {}
	OctreeNode(BoundingBox<T>& bBox) : weight(0), centerOfMass(Point<T>(bBox.getCenter().getDimensions())), children({}), bBox(bBox) {}

	inline bool isLeaf() const {
		return children.size() == 0;
	}

	inline bool isEmpty() const {
		return weight == 0;
	}

	inline bool contains(const Point<T>& point) const {
		return bBox.contains(point);
	}

	void computeCenterOfMass() {
		if (!isLeaf()) {
			centerOfMass.scale(1.0/(double) weight);

			// remove empty childs
			children.erase(std::remove_if(children.begin(), children.end(), [&](OctreeNode<T>& child){return child.isEmpty();}), children.end());

			for (auto &child : children) {
				child.computeCenterOfMass();
			}
		}
	}

	void compress() {
		if (!isLeaf()) {
			count numNonEmptyChilren = 0;
			while (numNonEmptyChilren <= 1) {
				numNonEmptyChilren = 0;
				OctreeNode<T> lastChild;
				for (auto &child : children) {
					if (!child.isEmpty()) {
						numNonEmptyChilren++;
						lastChild = child;
					}
				}

				if (numNonEmptyChilren == 1) { // compress
					children = lastChild.children;
					bBox = lastChild.bBox;
				}
			}

			for (auto &child : children) {
				child.compress();
			}
		}
	}

	std::string toString() {
		std::string str;
		str += bBox.getCenter().toString() + " sL=" + std::to_string(bBox.getSideLength()) + "w=" + std::to_string(weight);
		str += "(";
		for (auto &child : children) {
			str += "[ ";
			str += child.toString();
			str += " ]";
		}
		str += ")";
		return str;
	}

	void split(count dimensions, count numChildren) {
		children = std::vector<OctreeNode<T>>(numChildren, OctreeNode<T>(bBox));
		for (index i = 0; i < numChildren; ++i) { // 0-bit => center - halfSideLength, 1-bit => center + halfSideLength, least-significant bit is lowest dimension
			children[i].bBox.setSideLength(bBox.getHalfSideLength());
			for (index d = 0; d < dimensions; ++d) {
				if (((i & ~(~0 << (d+1))) >> d) == 0) { // 0-bit
					children[i].bBox.getCenter()[d] -= children[i].bBox.getHalfSideLength();
				} else {
					children[i].bBox.getCenter()[d] += children[i].bBox.getHalfSideLength();
				}
			}
		}
	}

	void addPoint(const Point<T>& point, count dimensions, count numChildrenPerNode) {
		if (weight == 0) { // empty leaf
			weight++; // we add a point
			centerOfMass = point;
		} else {
			if (isLeaf()) { // split the leaf!
				if (point.distance(centerOfMass) < 1e-3) {
					centerOfMass += point;
					weight++;
					return;
				}
				split(dimensions, numChildrenPerNode);
				for (auto &child : children) { // add leaf point to one of this node's new children
					if (child.contains(centerOfMass)) {
						child.addPoint(centerOfMass, dimensions, numChildrenPerNode);
						break;
					}
				}
			}

			for (auto &child : children) {
				if (child.contains(point)) {
					child.addPoint(point, dimensions, numChildrenPerNode);
					break;
				}
			}

			weight++;
			centerOfMass += point;
		}
	}
};


/**
 * @ingroup viz
 *
 * Implementation of a k-dimensional octree for the purpose of Barnes-Hut-Approximation.
 */
template<typename T>
class Octree {
public:
	Octree() = default;
	Octree(const std::vector<Vector>& points);

	void recomputeTree(const std::vector<Vector> &points);

	inline std::vector<std::pair<count, Point<T>>> approximateDistance(const Point<T>& p, const double theta) const {
		return approximateDistance(root, p, theta);
	}

	inline void approximateDistance(const Point<T>& p, const double theta, std::vector<std::pair<count, Point<T>>>& result) const {
		approximateDistance(root, p, theta, result);
	}

	template<typename L>
	inline void approximateDistance(const Point<T>& p, const double theta, L& handle) const {
		approximateDistance(root, p, theta*theta, handle);
	}

	std::string toString() {
		return root.toString();
	}

private:
	OctreeNode<T> root;
	count dimensions;
	count numChildrenPerNode;


	std::vector<std::pair<count, Point<T>>> approximateDistance(const OctreeNode<T>& node, const Point<T>& p, const double theta) const;
	void approximateDistance(const OctreeNode<T>& node, const Point<T>& p, const double theta, std::vector<std::pair<count, Point<T>>>& result) const;

	template<typename L>
	void approximateDistance(const OctreeNode<T>& node, const Point<T>& p, const double sqTheta, L& handle) const;
};

template<typename T>
Octree<T>::Octree(const std::vector<Vector>& points) {
	dimensions = points.size();
	Point<T> center(dimensions);
	numChildrenPerNode = pow(2, dimensions);
	T sideLength = 0;
	for (count d = 0; d < dimensions; ++d) {
		T minVal = points[d][0];
		T maxVal = points[d][0];

		for (index i = 1; i < points[d].getDimension(); ++i) {
			minVal = std::min(minVal, points[d][i]);
			maxVal = std::max(maxVal, points[d][i]);
		}

		sideLength = std::max(sideLength, fabs(maxVal - minVal) * 1.005); // add 0.5% to bounding box
		center[d] = (minVal + maxVal)/2.0;
	}

	root.bBox = {center, sideLength};

	for (index i = 0; i < points[0].getDimension(); ++i) {
		Point<T> p(points.size());
		for (count d = 0; d < dimensions; ++d) {
			p[d] = points[d][i];
		}

		root.addPoint(p, dimensions, numChildrenPerNode);
	}

	root.computeCenterOfMass();
//	root.compress();
}

template<typename T>
void Octree<T>::recomputeTree(const std::vector<Vector>& points) {
	root.children.clear();
	Point<T> center(dimensions);
	T sideLength = 0;
	for (count d = 0; d < dimensions; ++d) {
		T minVal = points[d][0];
		T maxVal = points[d][0];

		for (index i = 1; i < points[d].getDimension(); ++i) {
			minVal = std::min(minVal, points[d][i]);
			maxVal = std::max(maxVal, points[d][i]);
		}

		sideLength = std::max(sideLength, fabs(maxVal - minVal) * 1.005); // add 0.5% to bounding box
		center[d] = (minVal + maxVal)/2.0;
	}

	root.bBox = {center, sideLength};

	for (index i = 0; i < points[0].getDimension(); ++i) {
		Point<T> p(points.size());
		for (count d = 0; d < dimensions; ++d) {
			p[d] = points[d][i];
		}

		root.addPoint(p, dimensions, numChildrenPerNode);
	}

	root.computeCenterOfMass();
}

template<typename T>
std::vector<std::pair<count, Point<T>>> Octree<T>::approximateDistance(const OctreeNode<T>& node, const Point<T>& p, const double theta) const {
	if (node.isEmpty()) return {};
	if (node.isLeaf()) {
		if (node.centerOfMass == p) return {};
		return {std::make_pair(node.weight, node.centerOfMass)};
	} else {
		double dist = p.distance(node.centerOfMass);
		if (dist == 0 || node.bBox.getSideLength() <= theta * dist) {
			return {std::make_pair(node.weight, node.centerOfMass)};
		} else { // split further
			std::vector<std::pair<count, Point<T>>> points;
			points.reserve(node.weight);
			for (auto &child : node.children) {
				std::vector<std::pair<count, Point<T>>> childPoints = approximateDistance(child, p, theta);
				points.insert(points.end(), childPoints.begin(), childPoints.end());
			}
			return points;
		}
	}

	return {};
}

template<typename T>
void Octree<T>::approximateDistance(const OctreeNode<T>& node, const Point<T>& p, const double theta, std::vector<std::pair<count, Point<T>>>& result) const {
	if (node.isEmpty()) return;
	if (node.isLeaf()) {
		if (node.centerOfMass != p) {
			result.push_back(std::make_pair(node.weight, node.centerOfMass));
		}
	} else {
		double dist = p.distance(node.centerOfMass);
		if (dist == 0 || node.bBox.getSideLength() <= theta * dist) {
			result.push_back(std::make_pair(node.weight, node.centerOfMass));
		} else { // split further and recurse
			for (auto &child : node.children) {
				approximateDistance(child, p, theta, result);
			}
		}
	}
}

template<typename T> template<typename L>
void Octree<T>::approximateDistance(const OctreeNode<T>& node, const Point<T>& p, const double sqTheta, L& handle) const {
	//if (node.isEmpty()) return;
	if (!node.isLeaf()) {
		double sqDist = p.squaredDistance(node.centerOfMass);
		if (sqDist == 0 || node.bBox.getSqSideLength() <= sqTheta * sqDist) {
			handle(node.weight, node.centerOfMass, sqDist);
		} else { // split further and recurse
			for (auto &child : node.children) {
				approximateDistance(child, p, sqTheta, handle);
			}
		}
	} else if (node.centerOfMass != p) { // node is leaf and non-empty since octree only stores non-empty nodes
		handle(node.weight, node.centerOfMass, p.squaredDistance(node.centerOfMass));
	}
}


} /* namespace NetworKit */

#endif /* NETWORKIT_CPP_VIZ_OCTREE_H_ */

/*
 * Quadtree.h
 *
 *  Created on: 21.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#ifndef QUADTREECARTESIANEUCLID_H_
#define QUADTREECARTESIANEUCLID_H_

#include <vector>
#include <memory>
#include <cmath>
#include <omp.h>
#include <functional>
#include "QuadNodeCartesianEuclid.h"

namespace NetworKit {

template <class T>
class QuadtreeCartesianEuclid {
	friend class QuadTreeCartesianEuclidGTest;
public:
	/**
	 * @param maxR Radius of the managed area. Must be smaller than 1.
	 * @param theoreticalSplit If true, split cells to get the same area in each child cell. Default is false
	 * @param alpha dispersion Parameter of the point distribution. Only has an effect if theoretical split is true
	 * @param capacity How many points can inhabit a leaf cell before it is split up?
	 *
	 */
	QuadtreeCartesianEuclid(Point<double> lower = Point<double>({0.0, 0.0}), Point<double> upper = Point<double>({1.0, 1.0}), bool theoreticalSplit=false, count capacity=1000) {
		assert(lower.getDimensions() == upper.getDimensions());
		root = QuadNodeCartesianEuclid<T>(lower, upper, capacity, theoreticalSplit);
		this->lower = lower;
		this->upper = upper;
	}

	QuadtreeCartesianEuclid(const vector<Point<double> > &positions, const vector<T> &content, bool theoreticalSplit=false, count capacity=1000) {
		const count n = positions.size();
		assert(content.size() == n);
		assert(n > 0);

		this->dimension = positions[0].getDimensions();
		vector<double> lowerValue(dimension);
		vector<double> upperValue(dimension);
		for (index d = 0; d < dimension; d++) {
			lowerValue[d] = positions[0].at(d);
			upperValue[d] = positions[0].at(d);
		}

		for (Point<double> pos : positions) {
			assert(pos.getDimensions() == dimension);
			for (index d = 0; d < dimension; d++) {
				if (pos[d] < lowerValue[d]) lowerValue[d] = pos[d];
				if (pos[d] > upperValue[d]) upperValue[d] = pos[d];
			}
		}

		//the upper limit is open, so it needs to be above the points
		for (index d = 0; d < dimension; d++) {
			upperValue[d] = std::nextafter(upperValue[d], std::numeric_limits<double>::max());
		}
		this->lower = Point<double>(lowerValue);
		this->upper = Point<double>(upperValue);

		root = QuadNodeCartesianEuclid<T>(lower, upper, capacity, theoreticalSplit);
		for (index i = 0; i < n; i++) {
			assert(content[i] < n);
			root.addContent(content[i], positions[i]);
		}
	}

	/**
	 * @param newcomer content to be added at point x
	 * @param angle angular coordinate of x
	 * @param R radial coordinate of x
	 */
	void addContent(T newcomer, Point<double> pos) {
		root.addContent(newcomer, pos);
	}

	/**
	 * @param newcomer content to be removed at point x
	 * @param angle angular coordinate of x
	 * @param R radial coordinate of x
	 */
	bool removeContent(T toRemove, Point<double> pos) {
		return root.removeContent(toRemove, pos);
	}

	/**
	 * Get all elements, regardless of position
	 *
	 * @return vector<T> of elements
	 */
	vector<T> getElements() const {
		return root.getElements();
	}

	void extractCoordinates(vector<Point<double> > &posContainer) const {
		root.getCoordinates(posContainer);
	}

	void getElementsInEuclideanCircle(const Point<double> circleCenter, const double radius, vector<T> &circleDenizens) const {
		root.getElementsInEuclideanCircle(circleCenter, radius, circleDenizens);
	}

	template<typename L>
	count getElementsProbabilistically(Point<double> euQuery, L prob, vector<T> &circleDenizens) {
		return root.getElementsProbabilistically(euQuery, prob, circleDenizens);
	}

	void recount() {
		root.recount();
	}

	count size() const {
		return root.size();
	}

	count height() const {
		return root.height();
	}

	count countLeaves() const {
		return root.countLeaves();
	}

	index indexSubtree(index nextID) {
		return root.indexSubtree(nextID);
	}

	index getCellID(Point<double> pos) const {
		return root.getCellID(pos);
	}

	void reindex() {
		#pragma omp parallel
		{
			#pragma omp single nowait
			{
				root.reindex(0);
			}
		}
	}

	/**
	 * trims the vectors used to hold the content in the leaf cells. Reduces memory usage, makes changes slower
	 */
	void trim() {
		root.trim();
	}

private:
	QuadNodeCartesianEuclid<T> root;
	Point<double> lower;
	Point<double> upper;
	count dimension;
};
}

#endif /* QUADTREE_H_ */

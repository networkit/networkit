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
	QuadtreeCartesianEuclid() {
		root = QuadNodeCartesianEuclid<T>();
	}

	/**
	 * @param maxR Radius of the managed area. Must be smaller than 1.
	 * @param theoreticalSplit If true, split cells to get the same area in each child cell. Default is false
	 * @param alpha dispersion Parameter of the point distribution. Only has an effect if theoretical split is true
	 * @param capacity How many points can inhabit a leaf cell before it is split up?
	 *
	 */
	QuadtreeCartesianEuclid(Point2D<double> lowerLeft, Point2D<double> upperRight,bool theoreticalSplit=false, count capacity=1000) {
		root = QuadNodeCartesianEuclid<T>(lowerLeft, upperRight, capacity, theoreticalSplit);
		lower = lowerLeft;
		upper = upperRight;
	}

	QuadtreeCartesianEuclid(const vector<Point2D<double> > &positions, const vector<T> &content, bool theoreticalSplit=false, count capacity=1000) {
		const count n = positions.size();
		assert(content.size() == n);

		for (Point2D<double> pos : positions) {
			if (pos.getX() < lower.getX()) lower = Point2D<double>(pos.getX(), lower.getY());
			if (pos.getY() < lower.getY()) lower = Point2D<double>(lower.getX(), pos.getY());

			if (pos.getX() > upper.getX()) upper = Point2D<double>(pos.getX(), upper.getY());
			if (pos.getY() > upper.getY()) upper = Point2D<double>(upper.getX(), pos.getY());
		}

		//the upper limit is open, so it needs to be above the points
		upper = Point2D<double>(std::nextafter(upper.getX(), std::numeric_limits<double>::max()), std::nextafter(upper.getY(), std::numeric_limits<double>::max()));

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
	void addContent(T newcomer, Point2D<double> pos) {
		root.addContent(newcomer, pos);
	}

	/**
	 * @param newcomer content to be removed at point x
	 * @param angle angular coordinate of x
	 * @param R radial coordinate of x
	 */
	bool removeContent(T toRemove, Point2D<double> pos) {
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

	void extractCoordinates(vector<Point2D<double> > &posContainer) const {
		root.getCoordinates(posContainer);
	}

	void getElementsInEuclideanCircle(const Point2D<double> circleCenter, const double radius, vector<T> &circleDenizens) const {
		root.getElementsInEuclideanCircle(circleCenter, radius, false, circleDenizens);
	}

	count getElementsProbabilistically(Point2D<double> euQuery, std::function<double(double)> prob, vector<T> &circleDenizens) {
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

	index getCellID(Point2D<double> pos) const {
		return root.getCellID(pos);
	}


	void sortPointsInLeaves() {
		#pragma omp parallel
		{
			#pragma omp single nowait
			{
				root.sortPointsInLeaves();
			}
		}
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
	Point2D<double> lower;
	Point2D<double> upper;
};
}

#endif /* QUADTREE_H_ */

/*
 * QuadNode.h
 *
 *  Created on: 21.05.2014
 *      Author: moritz
 */

#ifndef QUADNODE_H_
#define QUADNODE_H_

#include <vector>
#include <algorithm>
#include <assert.h>
#include "../HyperbolicSpace.h"

using std::vector;

namespace NetworKit {

template <class T>
class QuadNode {
public:
	QuadNode() {
		leftAngle = 0;
		rightAngle = 2*M_PI;
		minR = 0;
		maxR = 10;//TODO: magic Number, careful.
		capacity = 20;
		isLeaf = true;
	}

	~QuadNode() {
		// TODO Auto-generated constructor stub

	}

	QuadNode(double leftAngle, double minR, double rightAngle, double maxR, unsigned capacity, double minDiameter) {
		this->leftAngle = leftAngle;
		this->minR = minR;
		this->rightAngle = rightAngle;
		this->maxR = maxR;
		this->capacity = 20;
		this->minRegion = minDiameter;
		isLeaf = true;
	}

	void addContent(T input, double angle, double R) {
		assert(this->responsible(angle, R));
		if (isLeaf) {
			if (content.size() + 1 < capacity ||  HyperbolicSpace::getDistance(leftAngle, minR, rightAngle, maxR) < minRegion) {
				content.push_back(input);
				angles.push_back(angle);
				radii.push_back(R);
			} else {
				//heavy lifting: split up!
				double middleAngle = (rightAngle - leftAngle) / 2 + leftAngle;
				double middleR = (maxR - minR) / 2 + minR;
				QuadNode southwest(leftAngle, minR, middleAngle, middleR, capacity, minRegion);
				QuadNode southeast(middleAngle, minR, rightAngle, middleR, capacity, minRegion);
				QuadNode northwest(leftAngle, middleR, middleAngle, maxR, capacity, minRegion);
				QuadNode northeast(middleAngle, middleR, rightAngle, maxR, capacity, minRegion);
				children = {southwest, southeast, northwest, northeast};

				isLeaf = false;
				for (uint i = 0; i < content.size(); i++) {
					this->addContent(content[i], angles[i], radii[i]);
				}
				content.clear();
				this->addContent(input, angle, R);
			}
		}
		else {
			assert(children.size() > 0);
			bool foundResponsibleChild = false;
			for (uint i = 0; i < children.size(); i++) {
				if (children[i].responsible(angle, R)) {
					assert(!foundResponsibleChild);//only one!
					children[i].addContent(input, angle, R);
					foundResponsibleChild = true;
				} else {
					//cout << "Not responsible for (" << angle << ", " << R << "). Borders are " << children[i].leftAngle << "-" << children[i].rightAngle << ", and " << children[i].minR << "-" << children[i].maxR << endl;
				}
			}
			assert(foundResponsibleChild);
		}
	}

	double distanceLowerBound(double angle, double R) {
		double nearestR = R;
		if (nearestR < minR) nearestR = minR;
		if (nearestR >= maxR) nearestR = maxR;
		if (angle >= this->leftAngle && angle < this->rightAngle) {
			return std::abs(R-nearestR);
		}
		double leftDistance = HyperbolicSpace::getDistance(angle, R, this->leftAngle, nearestR);
		double rightDistance = HyperbolicSpace::getDistance(angle, R, this->rightAngle, nearestR);
		return std::min(leftDistance, rightDistance);
	}

	double distanceUpperBound(double angle, double R) {
		double leftLower = HyperbolicSpace::getDistance(angle, R, this->leftAngle, this->minR);
		double rightLower = HyperbolicSpace::getDistance(angle, R, this->rightAngle, this->minR);
		double leftUpper = HyperbolicSpace::getDistance(angle, R, this->leftAngle, this->maxR);
		double rightUpper = HyperbolicSpace::getDistance(angle, R, this->rightAngle, this->maxR);
		return std::max(std::max(leftLower, rightLower), std::max(leftUpper, rightUpper));
	}

	bool responsible(double angle, double R) {
		return (angle >= leftAngle && angle < rightAngle && R >= minR && R < maxR);
	}

	std::vector<T> getElements() {
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

	std::vector<T> getCloseElements(double angle, double R, double maxDistance) {
		std::vector<T> result;
		if (isLeaf) {
			if (this->distanceLowerBound(angle, R) < maxDistance) {
				for (uint i = 0; i < content.size(); i++) {
					if (HyperbolicSpace::getDistance(angle, R, angles[i], radii[i]) < maxDistance) {
						result.push_back(content[i]);
					}
				}
				}
		} else {
			for (uint i = 0; i < children.size(); i++) {
				QuadNode * child = &children[i];
				if (child->distanceLowerBound(angle, R) < maxDistance) {
					vector<T> subresult = child->getCloseElements(angle, R, maxDistance);
					result.insert(result.end(), subresult.begin(), subresult.end());
				}
			}
		}
		return result;
	}

private:
	double leftAngle;
	double rightAngle;
	double minR;
	double maxR;
	unsigned capacity;
	double minRegion;//the minimal region a QuadNode should cover. If it is smaller, don't bother splitting up.
	std::vector<QuadNode> children;
	std::vector<T> content;
	std::vector<double> angles;
	std::vector<double> radii;
	bool isLeaf;
};
}

#endif /* QUADNODE_H_ */

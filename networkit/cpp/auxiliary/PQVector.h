/*
 * PQVector.h
 *
 *  Created on: 21.09.2018
 *      Author: Eugenio Angriman
 */

#ifndef PQVECTOR_H_
#define PQVECTOR_H_

#include <algorithm>
#include <vector>

#include "../auxiliary/Log.h"

namespace Aux {
/*
 * A simple prio queue for Kadabra that allows inspecting the ranking.
 */
class PQVector {
private:
	std::vector<std::pair<uint64_t, double>> elements;
	std::vector<uint64_t> position;
	uint64_t virtualSize;
	const uint64_t size;
	const uint64_t maxKey;

public:
	PQVector(const uint64_t size, const uint64_t maxKey);
	void insert(const uint64_t newElement, const double newValue);
	double getValue(const uint64_t i) const { return elements[i].second; }
	uint64_t getElement(const uint64_t i) const { return elements[i].first; }
	uint64_t getSize() const { return virtualSize; }
	void clear();
};

inline PQVector::PQVector(const uint64_t size, const uint64_t maxKey)
		: size(size), maxKey(maxKey) {
	if (maxKey < size) {
		throw std::runtime_error("maxKey must be bigger than the size.");
	}
	clear();
}

inline void PQVector::clear() {
	elements.resize(size);
	position.resize(maxKey);
	uint64_t i;
	for (i = 0; i < size; ++i) {
		elements[i] = std::make_pair(i, 0.);
		position[i] = i;
	}
	for (i = size; i < maxKey; ++i) {
		position[i] = i;
	}
	virtualSize = 0;
}

inline void PQVector::insert(const uint64_t newElement, const double newValue) {
	uint64_t ub =
			std::upper_bound(
					elements.begin(), elements.begin() + virtualSize, newValue,
					[&](const double val, const std::pair<uint64_t, double> pair) {
						return val > pair.second;
					}) -
			elements.begin();

	uint64_t oldPos;
	// We assume that if the same key is inserted again, its value will be
	// greater than before (i.e., oldPos >= ub).
	if (ub < size) {
		oldPos = std::min(position[newElement], size - 1);
		if (position[newElement] < size) {
			assert(elements[ub].second <= newValue);
		}
		if (virtualSize < size && oldPos >= virtualSize) {
			++virtualSize;
		}

		if (ub < oldPos) {
			std::rotate(elements.begin() + ub, elements.begin() + oldPos,
									elements.begin() + oldPos + 1);

			if (oldPos == size - 1) {
				++position[elements[ub].first];
			}

			elements[ub] = std::make_pair(newElement, newValue);
			for (auto it = elements.begin() + ub + 1;
					 it < elements.begin() + oldPos + 1; ++it) {
				++position[(*it).first];
			}
		} else {
			elements[ub].second = newValue;
		}
	}
}
} // namespace Aux

#endif

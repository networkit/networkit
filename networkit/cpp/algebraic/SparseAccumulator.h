/*
 * SparseAccumulator.h
 *
 *  Created on: 14.05.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef SPARSEACCUMULATOR_H_
#define SPARSEACCUMULATOR_H_

#include <vector>
#include <algorithm>
#include <cassert>

namespace NetworKit {

/**
 * The SparseAccumulator class represents the sparse accumulator datastructure as described in Kepner, Jeremy, and John Gilbert, eds.
 * Graph algorithms in the language of linear algebra. Vol. 22. SIAM, 2011. It is used as temporal storage for efficient computations on matrices.
 */
class SparseAccumulator {
private:
	/** row indicator used to avoid resetting the occupied vector */
	count row;

protected:
	/** dense vector of doubles which stores the computed values */
	std::vector<double> values; // w

	/** dense vector of integers which stores whether a valid value is stored in values at each position */
	std::vector<count> occupied; // b

	/** unordered list to store the position of valid values in values vector */
	std::vector<index> indices; // LS

public:
	/**
	 * Constructs the SparseAccumulator with size @a size.
	 * @param size The size of the SparseAccumulator.
	 */
	SparseAccumulator(count size) : row(1), values(size), occupied(size, 0) {}

	/**
	 * Stores @a value at @a pos. If a valid value is already stored at @a pos then @value is added to that.
	 * @param value The value to store or add at @a pos in values.
	 * @param pos The position in values.
	 */
	void scatter(double value, index pos) {
		assert(pos < values.size());

		if (occupied[pos] < row) {
			values[pos] = value;
			occupied[pos] = row;
			indices.push_back(pos);
		} else {
			values[pos] += value;
		}
	}

	/**
	 * Calls @a handle for each non zero value of the current row.
	 * @note handle signature: handle(index row, index column, double value)
	 * @return The number of non zero values in the current row.
	 *
	 */
	template<typename L> count gather(L handle) {
		count nonZeros = 0;
		std::sort(indices.begin(), indices.end());
		for (index idx : indices) {
			handle(row - 1, idx, values[idx]);
			nonZeros++;
		}

		return nonZeros;
	}

	/**
	 * Sets the SparseAccumulator to the next row which invalidates all currently stored data.
	 */
	void increaseRow() {
		row++;
		indices.clear();
	}
};

} /* namespace NetworKit */

#endif /* SPARSEACCUMULATOR_H_ */

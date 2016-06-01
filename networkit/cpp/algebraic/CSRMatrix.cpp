/*
 * CSRMatrix.cpp
 *
 *  Created on: May 6, 2015
 *      Author: Michael
 */

#include "CSRMatrix.h"

#include <cassert>
#include <atomic>
#include "omp.h"

namespace NetworKit {

CSRMatrix::CSRMatrix() : rowIdx(0), columnIdx(0), nonZeros(0), nRows(0), nCols(0), isSorted(true), zero(0.0) {
}

CSRMatrix::CSRMatrix(const count dimension, const double zero) : rowIdx(dimension+1), columnIdx(0), nonZeros(0), nRows(dimension), nCols(dimension), isSorted(true), zero(zero) {
}

CSRMatrix::CSRMatrix(const count nRows, const count nCols, const double zero) : rowIdx(nRows+1), columnIdx(0), nonZeros(0), nRows(nRows), nCols(nCols), isSorted(true), zero(zero) {
}

CSRMatrix::CSRMatrix(const count dimension, const std::vector<Triplet>& triplets, const double zero, bool isSorted) : CSRMatrix(dimension, dimension, triplets, zero, isSorted) {
}

CSRMatrix::CSRMatrix(const count nRows, const count nCols, const std::vector<Triplet>& triplets, const double zero, bool isSorted) : nRows(nRows), nCols(nCols), isSorted(isSorted), zero(zero) {
	count nnz = triplets.size();
	rowIdx = std::vector<index>(nRows + 1, 0);
	columnIdx = std::vector<index>(nnz);
	nonZeros = std::vector<double>(nnz);

	for (index i = 0; i < nnz; ++i) {
		rowIdx[triplets[i].row]++;
	}

	for (index i = 0, prefixSum = 0; i < nRows; ++i) {
		count nnzInRow = rowIdx[i];
		rowIdx[i] = prefixSum;
		prefixSum += nnzInRow;
	}
	rowIdx[nRows] = nnz;

	for (index i = 0; i < nnz; ++i) {
		index row = triplets[i].row;
		index dest = rowIdx[row];

		columnIdx[dest] = triplets[i].column;
		nonZeros[dest] = triplets[i].value;

		rowIdx[row]++;
	}

	for (index i = 0, firstIdxOfRow = 0; i <= nRows; ++i) {
		index newRow = rowIdx[i];
		rowIdx[i] = firstIdxOfRow;
		firstIdxOfRow = newRow;
	}
}

CSRMatrix::CSRMatrix(const count nRows, const count nCols, const std::vector<std::vector<index>> &columnIdx, const std::vector<std::vector<double>> &values,  const double zero, bool isSorted) : nRows(nRows), nCols(nCols), isSorted(isSorted), zero(zero) {
	 rowIdx = std::vector<index>(nRows + 1, 0);
	 count nnz = columnIdx[0].size();
	 for (index i = 1; i < columnIdx.size(); ++i) {
		 rowIdx[i] = rowIdx[i-1] + columnIdx[i-1].size();
		 nnz += columnIdx[i].size();
	 }
	 rowIdx[nRows] = nnz;

	 this->columnIdx = std::vector<index>(nnz);
	 this->nonZeros = std::vector<double>(nnz);

#pragma omp parallel for
	 for (index i = 0; i < nRows; ++i) {
		 for (index k = 0; k < columnIdx[i].size(); ++k) {
			 this->columnIdx[rowIdx[i] + k] = columnIdx[i][k];
			 nonZeros[rowIdx[i] + k] = values[i][k];
		 }
	 }
}

CSRMatrix::CSRMatrix(const count nRows, const count nCols, const std::vector<index>& rowIdx, const std::vector<index>& columnIdx, const std::vector<double>& nonZeros,  const double zero, bool isSorted) : rowIdx(rowIdx), columnIdx(columnIdx), nonZeros(nonZeros), nRows(nRows), nCols(nCols), isSorted(isSorted), zero(zero) {
}

count CSRMatrix::nnzInRow(const index i) const {
	assert(i >= 0 && i < nRows);
	return rowIdx[i+1] - rowIdx[i];
}

count CSRMatrix::nnz() const {
	return nonZeros.size();
}

double CSRMatrix::operator()(const index i, const index j) const {
	assert(i >= 0 && i < nRows);
	assert(j >= 0 && j < nCols);

	if (rowIdx[i] == rowIdx[i+1]) return zero; // no non-zero value is present in this row

	double value = zero;
//	if (!sorted()) {
		for (index k = rowIdx[i]; k < rowIdx[i+1]; ++k) {
			if (columnIdx[k] == j) {
				value = nonZeros[k];
				break;
			}
		}
//	} else {
//		index colIdx = binarySearchColumns(rowIdx[i], rowIdx[i+1]-1, j);
//		if (rowIdx[i] <= colIdx && colIdx < rowIdx[i+1] && columnIdx[colIdx] == j) {
//			value = nonZeros[colIdx];
//		}
//	}

	return value;
}

void CSRMatrix::setValue(const index i, const index j, const double value) {
	assert(i < nRows);
	assert(j < nCols);

	index colIdx = none;
	if (!sorted()) {
		for (index k = rowIdx[i]; k < rowIdx[i+1]; ++k) {
			if (columnIdx[k] == j) {
				colIdx = k;
			}
		}
	} else {
		colIdx = binarySearchColumns(rowIdx[i], rowIdx[i+1]-1, j);
	}

	if (colIdx != none && colIdx >= rowIdx[i] && columnIdx[colIdx] == j) { // the matrix already has an entry at (i,j) => replace it
		nonZeros[colIdx] = value;
	} else { // create a new non-zero entry at (i,j)
		if (!sorted()) {
			columnIdx.emplace(std::next(columnIdx.begin(), rowIdx[i+1]), j);
			nonZeros.emplace(std::next(nonZeros.begin(), rowIdx[i+1]), value);
		} else {
			if (colIdx < rowIdx[i]) { // emplace the value in front of all other values of row i
				columnIdx.emplace(std::next(columnIdx.begin(), rowIdx[i]), value);
				nonZeros.emplace(std::next(nonZeros.begin(), rowIdx[i]), value);
			} else {
				columnIdx.emplace(std::next(columnIdx.begin(), rowIdx[i] + colIdx+1), value);
				nonZeros.emplace(std::next(nonZeros.begin(), rowIdx[i] + colIdx+1), value);
			}
		}

		// update rowIdx
		for (index k = i+1; k < rowIdx.size(); ++k) {
			rowIdx[k]++;
		}
	}
}

void CSRMatrix::quicksort(index left, index right) {
	if (left >= right) return;
	index pivotIdx = partition(left, right);
	if (pivotIdx != 0) {
		quicksort(left, pivotIdx-1);
	}
	quicksort(pivotIdx+1, right);
}

index CSRMatrix::partition(index left, index right) {
	index mid = (left + right) / 2;
	index pivot = columnIdx[mid];
	std::swap(columnIdx[mid], columnIdx[right]);
	std::swap(nonZeros[mid], nonZeros[right]);

	index i = left;
	for (index j = left; j < right; ++j) {
		if (columnIdx[j] <= pivot) {
			std::swap(columnIdx[i], columnIdx[j]);
			std::swap(nonZeros[i], nonZeros[j]);
			++i;
		}
	}
	std::swap(columnIdx[i], columnIdx[right]);
	std::swap(nonZeros[i], nonZeros[right]);
	return i;
}

index CSRMatrix::binarySearchColumns(index left, index right, index j) const {
	assert(sorted());
	if (left > right) return right; // return the index immediately left to j if it would be present
	index mid = (left + right) / 2;
	if (columnIdx[mid] == j) {
		return mid;
	} else if (columnIdx[mid] > j) {
		return binarySearchColumns(left, mid-1, j);
	} else {
		return binarySearchColumns(mid+1, right, j);
	}
}

void CSRMatrix::sort() {
#pragma omp parallel for schedule(guided)
	for (index i = 0; i < nRows; ++i) {
		if (rowIdx[i+1] - rowIdx[i] > 1) {
			quicksort(rowIdx[i], rowIdx[i+1]-1);
		}

	}
	assert(sorted());

	isSorted = true;
}

bool CSRMatrix::sorted() const {
#ifndef NDEBUG
	bool sorted = true;
#pragma omp parallel for
	for (index i = 0; i < nRows; ++i) {
		for (index j = rowIdx[i]; j < rowIdx[i+1]; ++j) {
			for (index k = j+1; k < rowIdx[i+1]; ++k) {
				if (columnIdx[j] > columnIdx[k]) {
					sorted = false;
					break;
				}
			}
		}
	}

	return sorted;
#endif

	return isSorted;
}

Vector CSRMatrix::row(const index i) const {
	assert(i >= 0 && i < nRows);

	Vector row(numberOfColumns(), zero, true);
	parallelForNonZeroElementsInRow(i, [&](index j, double value) {
		row[j] = value;
	});

	return row;
}

Vector CSRMatrix::column(const index j) const {
	assert(j >= 0 && j < nCols);

	Vector column(numberOfRows(), getZero());
#pragma omp parallel for
	for (node i = 0; i < numberOfRows(); ++i) {
		column[i] = (*this)(i,j);
	}

	return column;
}

Vector CSRMatrix::diagonal() const {
	Vector diag(std::min(nRows, nCols), zero);

	if (sorted()) {
#pragma omp parallel for
		for (index i = 0; i < diag.getDimension(); ++i) {
			if (rowIdx[i] == rowIdx[i+1]) continue; // no entry in row i
			index left = rowIdx[i];
			index right = rowIdx[i+1]-1;
			index mid = (left + right) / 2;
			while (left <= right) {
				if (columnIdx[mid] == i) {
					diag[i] = nonZeros[mid];
					break;
				}

				if (columnIdx[mid] < i) {
					left = mid+1;
				} else {
					right = mid-1;
				}

				mid = (left + right) / 2;
			}
		}
	} else {
#pragma omp parallel for
		for (index i = 0; i < diag.getDimension(); ++i) {
			diag[i] = (*this)(i,i);
		}
	}

	return diag;
}

CSRMatrix CSRMatrix::operator+(const CSRMatrix &other) const {
	assert(nRows == other.nRows && nCols == other.nCols);
	return CSRMatrix::binaryOperator(*this, other, [](double val1, double val2) {return val1 + val2;});
}

CSRMatrix& CSRMatrix::operator+=(const CSRMatrix &other) {
	assert(nRows == other.nRows && nCols == other.nCols);
	*this = CSRMatrix::binaryOperator(*this, other, [](double val1, double val2) {return val1 + val2;});
	return *this;
}

CSRMatrix CSRMatrix::operator-(const CSRMatrix &other) const {
	assert(nRows == other.nRows && nCols == other.nCols);
	return CSRMatrix::binaryOperator(*this, other, [](double val1, double val2) {return val1 - val2;});
}

CSRMatrix& CSRMatrix::operator-=(const CSRMatrix &other) {
	assert(nRows == other.nRows && nCols == other.nCols);
	*this = CSRMatrix::binaryOperator(*this, other, [](double val1, double val2) {return val1 - val2;});
	return *this;
}

CSRMatrix CSRMatrix::operator*(const double &scalar) const {
	return CSRMatrix(*this) *= scalar;
}

CSRMatrix& CSRMatrix::operator*=(const double &scalar) {
#pragma omp parallel for
	for (index k = 0; k < nonZeros.size(); ++k) {
		nonZeros[k] *= scalar;
	}

	return *this;
}

Vector CSRMatrix::operator*(const Vector &vector) const {
	assert(!vector.isTransposed());
	assert(nCols == vector.getDimension());

	Vector result(nRows, zero);
#pragma omp parallel for
	for (index i = 0; i < numberOfRows(); ++i) {
		double sum = zero;
		for (index cIdx = rowIdx[i]; cIdx < rowIdx[i+1]; ++cIdx) {
			sum += nonZeros[cIdx] * vector[columnIdx[cIdx]];
		}
		result[i] = sum;
	}

	return result;
}

CSRMatrix CSRMatrix::operator*(const CSRMatrix &other) const {
	assert(nCols == other.nRows);

	std::vector<index> rowIdx(numberOfRows()+1, 0);
	std::vector<index> columnIdx;
	std::vector<double> nonZeros;

#pragma omp parallel
	{
		std::vector<int64_t> marker(other.numberOfColumns(), -1);
		count numThreads = omp_get_num_threads();
		index threadId = omp_get_thread_num();

		count chunkSize = (numberOfRows() + numThreads - 1) / numThreads;
		index chunkStart = threadId * chunkSize;
		index chunkEnd = std::min(numberOfRows(), chunkStart + chunkSize);

		for (index i = chunkStart; i < chunkEnd; ++i) {
			for (index jA = this->rowIdx[i]; jA < this->rowIdx[i+1]; ++jA) {
				index k = this->columnIdx[jA];
				for (index jB = other.rowIdx[k]; jB < other.rowIdx[k+1]; ++jB) {
					index j = other.columnIdx[jB];
					if (marker[j] != (int64_t) i) {
						marker[j] = i;
						++rowIdx[i+1];
					}
				}
			}
		}

		std::fill(marker.begin(), marker.end(), -1);

#pragma omp barrier
#pragma omp single
		{
			for (index i = 0; i < numberOfRows(); ++i) {
				rowIdx[i+1] += rowIdx[i];
			}

			columnIdx = std::vector<index>(rowIdx[numberOfRows()]);
			nonZeros = std::vector<double>(rowIdx[numberOfRows()]);
		}

		for (index i = chunkStart; i < chunkEnd; ++i) {
			index rowBegin = rowIdx[i];
			index rowEnd = rowBegin;

			for (index jA = this->rowIdx[i]; jA < this->rowIdx[i+1]; ++jA) {
				index k = this->columnIdx[jA];
				double valA = this->nonZeros[jA];

				for (index jB = other.rowIdx[k]; jB < other.rowIdx[k+1]; ++jB) {
					index j = other.columnIdx[jB];
					double valB = other.nonZeros[jB];

					if (marker[j] < (int64_t) rowBegin) {
						marker[j] = rowEnd;
						columnIdx[rowEnd] = j;
						nonZeros[rowEnd] = valA * valB;
						++rowEnd;
					} else {
						nonZeros[marker[j]] += valA * valB;
					}
				}
			}
		}
	}

	CSRMatrix result(numberOfRows(), other.numberOfColumns(), rowIdx, columnIdx, nonZeros);
	result.sort();
	return result;
}

CSRMatrix CSRMatrix::operator/(const double &divisor) const {
	return CSRMatrix(*this) /= divisor;
}

CSRMatrix& CSRMatrix::operator/=(const double &divisor) {
	return *this *= 1.0 / divisor;
}

CSRMatrix CSRMatrix::subMatrix(const std::vector<index> &rows, const std::vector<index> &columns) const {
	index invalid = std::numeric_limits<index>::max();
	std::vector<index> columnMapping(numberOfColumns(), invalid);
	std::vector<index> rowIdx(rows.size() + 1, 0);

#pragma omp parallel for
	for (index j = 0; j < columns.size(); ++j) {
		columnMapping[columns[j]] = j;
	}


#pragma omp parallel for
	for (index i = 0; i < rows.size(); ++i) {
		forNonZeroElementsInRow(rows[i], [&](index j, double val) {
			if (columnMapping[j] != invalid) {
				rowIdx[i+1]++;
			}
		});
	}

	for (index i = 0; i < rows.size(); ++i) {
		rowIdx[i+1] += rowIdx[i];
	}

	count nnz = rowIdx[rows.size()];
	std::vector<index> columnIdx(nnz);
	std::vector<double> nonZeros(nnz);

#pragma omp parallel for
	for (index i = 0; i < rows.size(); ++i) {
		index cIdx = rowIdx[i];
		forNonZeroElementsInRow(rows[i], [&](index j, double val) {
			if (columnMapping[j] != invalid) { // column is present in submatrix
				columnIdx[cIdx] = columnMapping[j];
				nonZeros[cIdx] = val;
				cIdx++;
			}
		});
	}

	return CSRMatrix(rows.size(), columns.size(), rowIdx, columnIdx, nonZeros, getZero(), sorted());
}

CSRMatrix CSRMatrix::mTmMultiply(const CSRMatrix &A, const CSRMatrix &B) {
	assert(A.nRows == B.nRows);

	std::vector<std::vector<index>> columnIdx(A.numberOfColumns());
	std::vector<std::vector<double>> values(A.numberOfColumns());

	for (index k = 0; k < A.numberOfRows(); ++k) {
		A.forNonZeroElementsInRow(k, [&](index i, double vA) {
			B.forNonZeroElementsInRow(k, [&](index j, double vB) {
				bool found = false;
				for (index l = 0; l < columnIdx[i].size(); ++l) {
					if (columnIdx[i][l] == j) {
						values[i][l] += vA * vB;
						found = true;
						break;
					}
				}

				if (!found) {
					columnIdx[i].push_back(j);
					values[i].push_back(vA * vB);
				}
			});
		});
	}

	return CSRMatrix(A.nCols, B.nCols, columnIdx, values);
}

CSRMatrix CSRMatrix::mmTMultiply(const CSRMatrix &A, const CSRMatrix &B) {
	assert(A.nCols == B.nCols);

	std::vector<std::vector<index>> columnIdx(A.numberOfRows());
	std::vector<std::vector<double>> values(A.numberOfRows());

	for (index i = 0; i < A.numberOfRows(); ++i) {
		A.forNonZeroElementsInRow(i, [&](index k, double vA) {
			for (index j = 0; j < B.numberOfRows(); ++j) {
				double vB = B(j,k);
				if (vB != A.zero) {
					bool found = false;
					for (index l = 0; l < columnIdx[i].size(); ++l) {
						if (columnIdx[i][l] == j) {
							values[i][l] += vA * vB;
							found = true;
							break;
						}
					}

					if (!found) {
						columnIdx[i].push_back(j);
						values[i].push_back(vA * vB);
					}
				}
			}
		});
	}

	return CSRMatrix(A.nRows, B.nRows, columnIdx, values);
}

Vector CSRMatrix::mTvMultiply(const CSRMatrix &matrix, const Vector &vector) {
	assert(matrix.nRows == vector.getDimension() && !vector.isTransposed());

	Vector result(matrix.numberOfColumns(), 0.0);
	for (index k = 0; k < matrix.numberOfRows(); ++k) {
		matrix.forNonZeroElementsInRow(k, [&](index j, double value) {
			result[j] += value * vector[k];
		});
	}

	return result;
}

CSRMatrix CSRMatrix::transpose() const {
	std::vector<index> rowIdx(numberOfColumns()+1);
	for (index i = 0; i < nnz(); ++i) {
		++rowIdx[columnIdx[i]+1];
	}

	for (index i = 0; i < numberOfColumns(); ++i) {
		rowIdx[i+1] += rowIdx[i];
	}

	std::vector<index> columnIdx(rowIdx[numberOfColumns()]);
	std::vector<double> nonZeros(rowIdx[numberOfColumns()]);

	for (index i = 0; i < numberOfRows(); ++i) {
		for (index j = this->rowIdx[i]; j < this->rowIdx[i+1]; ++j) {
			index colIdx = this->columnIdx[j];
			columnIdx[rowIdx[colIdx]] = i;
			nonZeros[rowIdx[colIdx]] = this->nonZeros[j];
			++rowIdx[colIdx];
		}
	}
	index shift = 0;
	for (index i = 0; i < numberOfColumns(); ++i) {
		index temp = rowIdx[i];
		rowIdx[i] = shift;
		shift = temp;
	}
	rowIdx[numberOfColumns()] = nonZeros.size();

	return CSRMatrix(nCols, nRows, rowIdx, columnIdx, nonZeros);
}

CSRMatrix CSRMatrix::adjacencyMatrix(const Graph &graph) {
	count nonZeros = graph.isDirected()? graph.numberOfEdges() : graph.numberOfEdges() * 2;
	std::vector<Triplet> triples(nonZeros);
	index index = 0;
	graph.forEdges([&](node i, node j, double val) {
		triples[index] = {i,j,val};
		index++;
		if (!graph.isDirected() && i != j) {
			triples[index] = {i,j,val};
			index++;
		}
	});

	return CSRMatrix(graph.upperNodeIdBound(), triples);
}

CSRMatrix CSRMatrix::diagonalMatrix(const Vector& diagonalElements) {
	count nRows = diagonalElements.getDimension();
	count nCols = diagonalElements.getDimension();
	std::vector<index> rowIdx(nRows+1, 0);
	std::iota(rowIdx.begin(), rowIdx.end(), 0);
	std::vector<index> columnIdx(nCols);
	std::vector<double> nonZeros(nCols);

#pragma omp parallel for
	for (index j = 0; j < nCols; ++j) {
		columnIdx[j] = j;
		nonZeros[j] = diagonalElements[j];
	}

	return CSRMatrix(nRows, nCols, rowIdx, columnIdx, nonZeros);
}

CSRMatrix CSRMatrix::incidenceMatrix(const Graph& graph) {
	if (!graph.hasEdgeIds()) throw std::runtime_error("Graph has no edge Ids. Index edges first by calling graph.indexEdges()");
	std::vector<Triplet> triples(2*graph.numberOfEdges()); // each edge produces two entries in the matrix

	index idx = 0;
	if (graph.isDirected()) {
		graph.forEdges([&](node u, node v, edgeweight weight, edgeid edgeId) {
			if (u != v) {
				edgeweight w = sqrt(weight);
				triples[idx++] = {u, edgeId, w};
				triples[idx++] = {v, edgeId, -w};
			}
		});
	} else {
		graph.forEdges([&](node u, node v, edgeweight weight, edgeid edgeId){
			if (u != v) {
				edgeweight w = sqrt(weight);
				if (u < v) { // orientation: small node number -> great node number
					triples[idx++] = {u, edgeId, w};
					triples[idx++] = {v, edgeId, -w};
				} else {
					triples[idx++] = {u, edgeId, -w};
					triples[idx++] = {v, edgeId, w};
				}
			}
		});
	}

	return CSRMatrix(graph.upperNodeIdBound(), graph.upperEdgeIdBound(), triples);
}

CSRMatrix CSRMatrix::laplacianMatrix(const Graph &graph) {
	std::vector<Triplet> triples;

	graph.forNodes([&](const index i){
		double weightedDegree = 0.0;
		graph.forNeighborsOf(i, [&](const index j, double weight) { // - adjacency matrix
			if (i != j) { // exclude diagonal since this would be subtracted by the adjacency weight
				weightedDegree += weight;
			}

			triples.push_back({i,j,-weight});
		});

		triples.push_back({i,i, weightedDegree}); // degree matrix
	});

	return CSRMatrix(graph.upperNodeIdBound(), triples);
}

CSRMatrix CSRMatrix::normalizedLaplacianMatrix(const Graph& graph) {
	std::vector<Triplet> triples;

	std::vector<double> weightedDegrees(graph.upperNodeIdBound(), 0.0);
	graph.parallelForNodes([&](const node u) {
		weightedDegrees[u] = graph.weightedDegree(u);
	});

	graph.forNodes([&](const node i){
		graph.forNeighborsOf(i, [&](const node j, double weight){
			if (i != j) {
				triples.push_back({i, j, -weight/sqrt(weightedDegrees[i] * weightedDegrees[j])});
			}
		});

		if (weightedDegrees[i] != 0.0) {
			if (graph.isWeighted()) {
				triples.push_back({i, i, 1-(graph.weight(i, i)) / weightedDegrees[i]});
			} else {
				triples.push_back({i, i, 1});
			}
		}
	});

	return CSRMatrix(graph.upperNodeIdBound(), triples);
}



} /* namespace NetworKit */

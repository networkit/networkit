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

/** Floating point epsilon to use in comparisons. */
constexpr double EPSILON = 1e-9;

CSRMatrix::CSRMatrix() : rowIdx(0), columnIdx(0), nonZeros(0), nRows(0), nCols(0), isSorted(true) {
}

CSRMatrix::CSRMatrix(const count nRows, const count nCols, const std::vector<std::pair<index, index>> &positions, const std::vector<double> &values, bool isSorted) : nRows(nRows), nCols(nCols), isSorted(isSorted) {
	count nnz = values.size();
	rowIdx = std::vector<index>(nRows + 1, 0);
	columnIdx = std::vector<index>(nnz);
	nonZeros = std::vector<double>(nnz);

	for (index i = 0; i < nnz; ++i) {
		rowIdx[positions[i].first]++;
	}

	for (index i = 0, prefixSum = 0; i < nRows; ++i) {
		count nnzInRow = rowIdx[i];
		rowIdx[i] = prefixSum;
		prefixSum += nnzInRow;
	}
	rowIdx[nRows] = nnz;

	for (index i = 0; i < nnz; ++i) {
		index row = positions[i].first;
		index dest = rowIdx[row];

		columnIdx[dest] = positions[i].second;
		nonZeros[dest] = values[i];

		rowIdx[row]++;
	}

	for (index i = 0, firstIdxOfRow = 0; i <= nRows; ++i) {
		index newRow = rowIdx[i];
		rowIdx[i] = firstIdxOfRow;
		firstIdxOfRow = newRow;
	}
}

CSRMatrix::CSRMatrix(const count nRows, const count nCols, const std::vector<Triple> &triples, bool isSorted) : nRows(nRows), nCols(nCols), isSorted(isSorted) {
	count nnz = triples.size();
	rowIdx = std::vector<index>(nRows + 1, 0);
	columnIdx = std::vector<index>(nnz);
	nonZeros = std::vector<double>(nnz);

	for (index i = 0; i < nnz; ++i) {
		rowIdx[triples[i].row]++;
	}

	for (index i = 0, prefixSum = 0; i < nRows; ++i) {
		count nnzInRow = rowIdx[i];
		rowIdx[i] = prefixSum;
		prefixSum += nnzInRow;
	}
	rowIdx[nRows] = nnz;

	for (index i = 0; i < nnz; ++i) {
		index row = triples[i].row;
		index dest = rowIdx[row];

		columnIdx[dest] = triples[i].column;
		nonZeros[dest] = triples[i].value;

		rowIdx[row]++;
	}

	for (index i = 0, firstIdxOfRow = 0; i <= nRows; ++i) {
		index newRow = rowIdx[i];
		rowIdx[i] = firstIdxOfRow;
		firstIdxOfRow = newRow;
	}
}

CSRMatrix::CSRMatrix(const count nRows, const count nCols, const std::vector<std::vector<index>> &columnIdx, const std::vector<std::vector<double>> &values, bool isSorted) : nRows(nRows), nCols(nCols), isSorted(isSorted) {
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

CSRMatrix::CSRMatrix(const count nRows, const count nCols, const std::vector<index> &rowIdx, const std::vector<index> &columnIdx, const std::vector<double> &nonZeros, bool isSorted) : rowIdx(rowIdx), columnIdx(columnIdx), nonZeros(nonZeros), nRows(nRows), nCols(nCols), isSorted(isSorted) {
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

	double value = 0.0;
	for (index k = rowIdx[i]; k < rowIdx[i+1]; ++k) {
		if (columnIdx[k] == j) {
			value = nonZeros[k];
			break;
		}
	}

	return value;
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

void CSRMatrix::sort() {
#pragma omp parallel for schedule(guided)
	for (index i = 0; i < nRows; ++i) {
		if (rowIdx[i+1] - rowIdx[i] > 1) {
			quicksort(rowIdx[i], rowIdx[i+1]-1);
		}
	}

	isSorted = true;
}

bool CSRMatrix::sorted() const {
#ifndef NDEBUG
	bool sorted = true;
#pragma omp parallel for
	for (index i = 0; i < nRows; ++i) {
		for (index j = rowIdx[i]+1; j < rowIdx[i+1]; ++j) {
			if (columnIdx[j-1] > columnIdx[j]) {
				sorted = false;
				break;
			}
		}
	}

	return sorted;
#endif

	return isSorted;
}

Vector CSRMatrix::row(const index i) const {
	assert(i >= 0 && i < nRows);

	Vector row(numberOfColumns(), 0.0, true);
	parallelForNonZeroElementsInRow(i, [&](index j, double value) {
		row[j] = value;
	});

	return row;
}

Vector CSRMatrix::column(const index j) const {
	assert(j >= 0 && j < nCols);

	Vector column(numberOfRows());
#pragma omp parallel for
	for (node i = 0; i < numberOfRows(); ++i) {
		column[i] = (*this)(i,j);
	}

	return column;
}

Vector CSRMatrix::diagonal() const {
	Vector diag(std::min(nRows, nCols), 0.0);

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

	Vector result(nRows, 0.0);
#pragma omp parallel for
	for (index i = 0; i < numberOfRows(); ++i) {
		double sum = 0.0;
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

//	std::vector<Triple> triples;
//
//	SparseAccumulator spa(numberOfRows());
//	for (index i = 0; i < numberOfRows(); ++i) {
//		forNonZeroElementsInRow(i, [&](index k, double val1) {
//			other.forNonZeroElementsInRow(k, [&](index j, double val2) {
//				spa.scatter(val1 * val2, j);
//			});
//		});
//
//		spa.gather([&](index i, index j, double value){
//			triples.push_back({i,j,value});
//		});
//
//		spa.increaseRow();
//	}
//
//	return CSRMatrix(nRows, other.nCols, triples, true);

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

	return CSRMatrix(rows.size(), columns.size(), rowIdx, columnIdx, nonZeros, sorted());
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
				if (vB != 0.0) {
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

CSRMatrix CSRMatrix::graphLaplacian(const Graph &graph) {
	std::vector<std::pair<index, index>> positions;
	std::vector<double> values;

	graph.forNodes([&](const index i){
		double weightedDegree = 0.0;

		double selfLoopWeight = 0.0;
		graph.forNeighborsOf(i, [&](const index j, double weight) { // - adjacency matrix
			if (j == i) {
				selfLoopWeight = weight;
			} else {
				positions.push_back(std::make_pair(i,j));
				values.push_back(-weight);
			}

			weightedDegree += weight;
		});

		positions.push_back(std::make_pair(i,i));
		values.push_back(weightedDegree - selfLoopWeight); // degree matrix
	});

	return CSRMatrix(graph.upperNodeIdBound(), graph.upperNodeIdBound(), positions, values);
}

CSRMatrix CSRMatrix::adjacencyMatrix(const Graph &graph) {
	int nonZeros = graph.isDirected()? graph.numberOfEdges() : graph.numberOfEdges() * 2;

	std::vector<std::pair<index, index>> positions(nonZeros);
	std::vector<double> values(nonZeros);

	int index = 0;
	graph.forEdges([&](node i, node j, double val) {
		positions[index] = std::make_pair(i,j);
		values[index] = val;
		index++;
		if (!graph.isDirected() && i != j) {
			positions[index] = std::make_pair(j,i);
			values[index] = val;
			index++;
		}
	});

	return CSRMatrix(graph.numberOfNodes(), graph.numberOfNodes(), positions, values);
}

Graph CSRMatrix::laplacianToGraph(const CSRMatrix &laplacian) {
	assert(isLaplacian(laplacian));
	Graph G(std::max(laplacian.numberOfRows(), laplacian.numberOfColumns()), true, false);
	laplacian.forNonZeroElementsInRowOrder([&](node u, node v, edgeweight weight) {
		if (u != v) { // exclude diagonal
			if (u < v) {
				G.addEdge(u, v, -weight);
			}
		}
	});

	return G;
}

Graph CSRMatrix::matrixToGraph(const CSRMatrix &matrix) {
	bool directed = !isSymmetric(matrix);
	Graph G(std::max(matrix.numberOfRows(), matrix.numberOfColumns()), true, directed);
	matrix.forNonZeroElementsInRowOrder([&](node u, node v, edgeweight weight) {
		if (directed || u <= v) {
			G.addEdge(u, v, weight);
		}
	});

	return G;
}

bool CSRMatrix::isSymmetric(const CSRMatrix &matrix) {
	bool output = true;
	matrix.forNonZeroElementsInRowOrder([&] (index i, index j, edgeweight w) {
		if (abs(matrix(j, i)-w) > EPSILON) {
			output = false;
		}
	});
	if (!output) INFO("not symmetric!");
	return output;
}

bool CSRMatrix::isSDD(const CSRMatrix &matrix) {
	if (!isSymmetric(matrix)) {
		return false;
	}

	/* Criterion: a_ii >= \sum_{j != i} a_ij */
	std::vector<double> row_sum(matrix.numberOfRows());
	matrix.parallelForNonZeroElementsInRowOrder([&] (node i, node j, double value) {
		if (i == j) {
			row_sum[i] += value;
		} else {
			row_sum[i] -= abs(value);
		}
	});

	return std::all_of(row_sum.begin(), row_sum.end(), [] (double val) {return val > -EPSILON;});
}

bool CSRMatrix::isLaplacian(const CSRMatrix &matrix) {
	if (!isSymmetric(matrix)) {
		return false;
	}

	/* Criterion: \forall_i \sum_j A_ij = 0  */
	std::vector<double> row_sum(matrix.numberOfRows());
	std::atomic<bool> right_sign(true);
	matrix.parallelForNonZeroElementsInRowOrder([&] (node i, node j, double value) {
		if (i != j && value > EPSILON) {
			right_sign = false;
		}
		row_sum[i] += value;
	});

	return right_sign && std::all_of(row_sum.begin(), row_sum.end(), [] (double val) {return abs(val) < EPSILON;});
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



} /* namespace NetworKit */

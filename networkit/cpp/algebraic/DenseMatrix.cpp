/*
 * DenseMatrix.cpp
 *
 *  Created on: Nov 25, 2015
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include <algorithm>

#include <networkit/algebraic/DenseMatrix.hpp>

namespace NetworKit {

DenseMatrix::DenseMatrix() : nRows(0), nCols(0), entries(std::vector<double>(0)), zero(0.0) {}

DenseMatrix::DenseMatrix(const count dimension, const double zero)
    : nRows(dimension), nCols(dimension), entries(dimension * dimension, zero), zero(zero) {}

DenseMatrix::DenseMatrix(const count nRows, const count nCols, const double zero)
    : nRows(nRows), nCols(nCols), entries(nRows * nCols, zero), zero(zero) {}

DenseMatrix::DenseMatrix(const count dimension, const std::vector<Triplet> &triplets,
                         const double zero)
    : DenseMatrix(dimension, dimension, triplets, zero) {}

DenseMatrix::DenseMatrix(const count nRows, const count nCols, const std::vector<Triplet> &triplets,
                         const double zero)
    : nRows(nRows), nCols(nCols), entries(nRows * nCols, zero), zero(zero) {
#pragma omp parallel for
    for (omp_index k = 0; k < static_cast<omp_index>(triplets.size()); ++k) {
        entries[triplets[k].row * nCols + triplets[k].column] = triplets[k].value;
    }
}

DenseMatrix::DenseMatrix(const count nRows, const count nCols, const std::vector<double> &entries,
                         const double zero)
    : nRows(nRows), nCols(nCols), entries(entries), zero(zero) {
    assert(entries.size() == nRows * nCols);
}

count DenseMatrix::nnzInRow(const index i) const {
    count nnz = 0;
    for (index offset = i * numberOfColumns(); offset < (i + 1) * numberOfColumns(); ++offset) {
        if (std::fabs(entries[offset] - zero) > FLOAT_EPSILON)
            nnz++;
    }
    return nnz;
}

count DenseMatrix::nnz() const {
    count nnz = 0;
    for (index k = 0; k < entries.size(); ++k) {
        if (std::fabs(entries[k] - zero) > FLOAT_EPSILON)
            nnz++;
    }
    return nnz;
}

double DenseMatrix::operator()(const index i, const index j) const {
    return entries[i * numberOfColumns() + j];
}

void DenseMatrix::setValue(const index i, const index j, const double value) {
    entries[i * numberOfColumns() + j] = value;
}

double &DenseMatrix::operator()(const index i, const index j) {
    return entries[i * numberOfColumns() + j];
}

Vector DenseMatrix::row(const index i) const {
    Vector row(numberOfColumns(), zero, true);
    index offset = i * numberOfColumns();
#pragma omp parallel for
    for (omp_index j = 0; j < static_cast<omp_index>(numberOfColumns()); ++j) {
        row[j] = entries[offset + j];
    }

    return row;
}

Vector DenseMatrix::column(const index j) const {
    Vector column(numberOfRows(), zero);
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(numberOfRows()); ++i) {
        column[i] = entries[i * numberOfColumns() + j];
    }

    return column;
}

Vector DenseMatrix::diagonal() const {
    Vector diagonal(std::min(numberOfRows(), numberOfColumns()), zero);
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(diagonal.getDimension()); ++i) {
        diagonal[i] = (*this)(i, i);
    }

    return diagonal;
}

DenseMatrix DenseMatrix::operator+(const DenseMatrix &other) const {
    assert(numberOfRows() == other.numberOfRows() && numberOfColumns() == other.numberOfColumns());
    return DenseMatrix::binaryOperator(*this, other,
                                       [](double val1, double val2) { return val1 + val2; });
}

DenseMatrix &DenseMatrix::operator+=(const DenseMatrix &other) {
    assert(numberOfRows() == other.numberOfRows() && numberOfColumns() == other.numberOfColumns());
    *this = DenseMatrix::binaryOperator(*this, other,
                                        [](double val1, double val2) { return val1 + val2; });
    return *this;
}

DenseMatrix DenseMatrix::operator-(const DenseMatrix &other) const {
    assert(numberOfRows() == other.numberOfRows() && numberOfColumns() == other.numberOfColumns());
    return DenseMatrix::binaryOperator(*this, other,
                                       [](double val1, double val2) { return val1 - val2; });
}

DenseMatrix &DenseMatrix::operator-=(const DenseMatrix &other) {
    assert(numberOfRows() == other.numberOfRows() && numberOfColumns() == other.numberOfColumns());
    *this = DenseMatrix::binaryOperator(*this, other,
                                        [](double val1, double val2) { return val1 - val2; });
    return *this;
}

DenseMatrix DenseMatrix::operator*(double scalar) const {
    return DenseMatrix(*this) *= scalar;
}

DenseMatrix &DenseMatrix::operator*=(double scalar) {
#pragma omp parallel for
    for (omp_index k = 0; k < static_cast<omp_index>(entries.size()); ++k) {
        entries[k] *= scalar;
    }

    return *this;
}

Vector DenseMatrix::operator*(const Vector &vector) const {
    assert(!vector.isTransposed());
    assert(numberOfColumns() == vector.getDimension());

    Vector result(numberOfRows(), zero);
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(numberOfRows()); ++i) {
        index offset = i * numberOfColumns();
        for (index k = offset, j = 0; k < offset + numberOfColumns(); ++k, ++j) {
            result[i] += entries[k] * vector[j];
        }
    }

    return result;
}

DenseMatrix DenseMatrix::operator*(const DenseMatrix &other) const {
    assert(numberOfColumns() == other.numberOfRows());
    std::vector<double> resultEntries(numberOfRows() * other.numberOfColumns());

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(numberOfRows()); ++i) {
        index offset = i * other.numberOfColumns();
        for (index k = 0; k < numberOfColumns(); ++k) {
            double val_i_k = (*this)(i, k);
            for (index j = 0; j < other.numberOfColumns(); ++j) {
                resultEntries[offset + j] += val_i_k * other(k, j);
            }
        }
    }

    return DenseMatrix(numberOfRows(), other.numberOfColumns(), resultEntries);
}

DenseMatrix DenseMatrix::operator/(double divisor) const {
    return DenseMatrix(*this) /= divisor;
}

DenseMatrix &DenseMatrix::operator/=(double divisor) {
    return *this *= 1.0 / divisor;
}

DenseMatrix DenseMatrix::transpose() const {
    DenseMatrix transposedMatrix(
        numberOfColumns(), numberOfRows(),
        std::vector<double>(numberOfRows() * numberOfColumns(), getZero()));
    forElementsInRowOrder(
        [&](index i, index j, double value) { transposedMatrix.setValue(j, i, value); });

    return transposedMatrix;
}

DenseMatrix DenseMatrix::extract(const std::vector<index> &rowIndices,
                                 const std::vector<index> &columnIndices) const {
    DenseMatrix result(rowIndices.size(), columnIndices.size(),
                       std::vector<double>(rowIndices.size() * columnIndices.size(), getZero()));
    for (index i = 0; i < rowIndices.size(); ++i) {
        for (index j = 0; j < columnIndices.size(); ++j) {
            double value = (*this)(rowIndices[i], columnIndices[j]);
            if (std::fabs(value - getZero()) > FLOAT_EPSILON) {
                result.setValue(i, j, value);
            }
        }
    }

    return result;
}

void DenseMatrix::assign(const std::vector<index> &rowIndices,
                         const std::vector<index> &columnIndices, const DenseMatrix &source) {
    assert(rowIndices.size() == source.numberOfRows());
    assert(columnIndices.size() == source.numberOfColumns());

    for (index i = 0; i < rowIndices.size(); ++i) {
        source.forElementsInRow(
            i, [&](index j, double value) { setValue(rowIndices[i], columnIndices[j], value); });
    }
}

DenseMatrix DenseMatrix::adjacencyMatrix(const Graph &graph, double zero) {
    DenseMatrix A(graph.upperNodeIdBound(), zero);
    graph.forEdges([&](node u, node v, edgeweight w) {
        A.setValue(u, v, w);
        if (!graph.isDirected()) { // add symmetric value at (v, u)
            A.setValue(v, u, w);
        }
    });

    return A;
}

DenseMatrix DenseMatrix::diagonalMatrix(const Vector &diagonalElements, double zero) {
    DenseMatrix D(diagonalElements.getDimension(), zero);
    for (index i = 0; i < diagonalElements.getDimension(); ++i) {
        D.setValue(i, i, diagonalElements[i]);
    }

    return D;
}

DenseMatrix DenseMatrix::incidenceMatrix(const Graph &graph, double zero) {
    if (!graph.hasEdgeIds())
        throw std::runtime_error(
            "Graph has no edge Ids. Index edges first by calling graph.indexEdges()");
    DenseMatrix I(graph.upperNodeIdBound(), graph.upperEdgeIdBound(), zero);
    if (graph.isDirected()) {
        graph.forEdges([&](node u, node v, edgeweight weight, edgeid edgeId) {
            if (u != v) {
                edgeweight w = std::sqrt(weight);
                I.setValue(u, edgeId, w);
                I.setValue(v, edgeId, -w);
            }
        });
    } else {
        graph.forEdges([&](node u, node v, edgeweight weight, edgeid edgeId) {
            if (u != v) {
                edgeweight w = std::sqrt(weight);
                if (u < v) { // orientation: small node number -> great node number
                    I.setValue(u, edgeId, w);
                    I.setValue(v, edgeId, -w);
                } else {
                    I.setValue(u, edgeId, -w);
                    I.setValue(v, edgeId, w);
                }
            }
        });
    }

    return I;
}

DenseMatrix DenseMatrix::laplacianMatrix(const Graph &graph, double zero) {
    DenseMatrix L(graph.numberOfNodes(), zero);
    graph.forNodes([&](const index i) {
        double weightedDegree = 0.0;

        graph.forNeighborsOf(i, [&](const index j, double weight) { // - adjacency matrix
            L.setValue(i, j, -weight);
            if (i != j) { // exclude weight of diagonal since it would be subtracted later
                weightedDegree += weight;
            }
        });

        L.setValue(i, i, weightedDegree); // degree matrix
    });

    return L;
}

void DenseMatrix::LUDecomposition(DenseMatrix &matrix) {
    assert(matrix.numberOfRows() == matrix.numberOfColumns());
    for (index k = 0; k < matrix.numberOfRows() - 1; ++k) {
        for (index i = k + 1; i < matrix.numberOfRows(); ++i) {
            matrix.setValue(i, k, matrix(i, k) / matrix(k, k));
            for (index j = k + 1; j < matrix.numberOfRows(); ++j) {
                matrix.setValue(i, j, matrix(i, j) - (matrix(i, k) * matrix(k, j)));
            }
        }
    }
}

Vector DenseMatrix::LUSolve(const DenseMatrix &LU, const Vector &b) {
    Vector x = b;

    for (index i = 0; i < LU.numberOfRows() - 1; ++i) { // forward substitution
        for (index j = i + 1; j < LU.numberOfRows(); ++j) {
            x[j] -= x[i] * LU(j, i);
        }
    }

    for (index i = LU.numberOfRows(); i-- > 0;) { // backward substitution
        x[i] /= LU(i, i);
        for (index j = 0; j < i; ++j) {
            x[j] -= x[i] * LU(j, i);
        }
    }

    return x;
}

} /* namespace NetworKit */

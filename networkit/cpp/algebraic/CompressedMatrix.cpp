/*
 * CompressedMatrix.cpp
 *
 *  Created on: 28.06.2014
 *      Author: dhoske
 */

#include "CompressedMatrix.h"
#include "../auxiliary/Parallelism.h"
#include <utility>
#include <omp.h>
#include <climits>

namespace NetworKit {

CompressedMatrix::CompressedMatrix(const Matrix& A)
    : nrows(A.numberOfRows()), row_idx(A.numberOfRows() + 1) {
  // Convert to CSR by visiting all of the rows
  std::vector<std::vector<std::pair<index, edgeweight>>> rows(nrows);
  A.parallelForNonZeroElementsInRowOrder([&] (index i, index j, edgeweight value) {
    rows[i].emplace_back(j, value);
  });

  // Sort rows for better cache efficiency
  for (index i = 0; i < nrows; ++i) {
    sort(rows[i].begin(), rows[i].end());
  }

  // Now construct CSR
  index idx = 0;
  row_idx[0] = 0;
  for (index i = 0; i < nrows; ++i) {
    for (const auto& entry: rows[i]) {
      cols.emplace_back(entry.first);
      entries.emplace_back(entry.second);
      idx++;
    }
    row_idx[i + 1] = idx;
  };
}

uint64_t CompressedMatrix::numberOfRows() const {
  return row_idx.size() - 1;
}

uint64_t CompressedMatrix::numberOfColumns() const {
  return numberOfRows();
}

double CompressedMatrix::operator()(index i, index j) const {
  assert(i < numberOfRows());
  for (index idx = row_idx[i]; idx < row_idx[i + 1]; ++idx) {
    if (cols[idx] == j) {
      return entries[idx];
    }
  }
  return 0.0;
}

Vector CompressedMatrix::operator*(const Vector& b) const {
  assert(row_idx.size() == b.getDimension() + 1);

  // Parallel Vector-matrix product
  Vector out(b.getDimension(), 0.0);
  #pragma omp parallel for schedule(dynamic, 2048)
  for (index row = 0; row < row_idx.size() - 1; ++row) {
    double cur_val = 0.0;
    for (index i = row_idx[row]; i < row_idx[row+1]; ++i) {
      cur_val += entries[i] * b[cols[i]];
    }
    out[row] = cur_val;
  }

  return out;
}

} /* namespace NetworKit */


/*
 * CompressedMatrix.h
 *
 *  Created on: 28.06.2014
 *      Author: dhoske
 */

#ifndef COMPRESSEDMATRIX_H_
#define COMPRESSEDMATRIX_H_

#include <cstdint>
#include <vector>

#include "Matrix.h"
#include "Vector.h"

namespace NetworKit {

/**
 * @ingroup algebraic
 * Matrix implementation using compressed sparse row storage. This matrix
 * implementation does not adhere to the complete matrix interface
 * (@see Matrix) and is quite incomplete. We basically only support
 * matrix-vector-products and access to entries.
 */
class CompressedMatrix {
public:
  /**
   * Constructs the compressed representation of @a A.
   */
  explicit CompressedMatrix(const Matrix& A);

  /**
   * Matrix-vector-product \f$Ax\f$. (LHS should explictly be
   * a CompressedMatrix, i.e. this is deliberately a member function.)
   */
  Vector operator*(const Vector& x) const;

  /**
   * Returns the number of rows of this matrix.
   */
  std::uint64_t numberOfRows() const;

  /**
   * Returns the number of columns of this matrix.
   */
  std::uint64_t numberOfColumns() const;

  /**
   * Returns the entry at position @a (i, j) of this matrix.
   */
  double operator()(index i, index j) const;

private:
  count nrows;
  std::vector<index>  row_idx;
  std::vector<double> entries;
  std::vector<index>  cols;
};

} /* namespace NetworKit */

#endif /* ADJACENCYMATRIX_H_ */

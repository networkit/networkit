/*
 * GraphReader.h
 *
 *  Created on: 17.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef MATRIXREADER_H_
#define MATRIXREADER_H_

#include "../algebraic/CSRMatrix.h"

namespace NetworKit {

/**
 * @ingroup io
 * Abstract base class for matrix readers.
 */
class MatrixReader {
public:
  virtual ~MatrixReader() = default;

  /**
   * Reads the matrix in @a path.
   */
  virtual CSRMatrix read(const std::string& path) = 0;

  /** only to be used by cython - this eliminates an unnecessary copy */
  CSRMatrix* _read(const std::string& path) {
    return new CSRMatrix{read(path)};
  };
};

}
#endif

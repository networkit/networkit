/*
 * MatrixMarketReader.h
 *
 *  Created on: 25.07.2014
 *      Author: dhoske
 */

#ifndef MATRIXMARKETREADER_H_
#define MATRIXMARKETREADER_H_

#include "MatrixReader.h"

namespace NetworKit {

/**
 * @ingroup io
 * Reader for the matrix market file format as described in
 * http://math.nist.gov/MatrixMarket/reports/MMformat.ps.
 *
 * Only allows real, symmetric matrices.
 */
class MatrixMarketReader: public NetworKit::MatrixReader {
public:

  MatrixMarketReader() = default;

  virtual CSRMatrix read(const std::string& path) override;

  /** Reads the matrix in @a in. */
  CSRMatrix read(std::istream& in);
};

}
#endif

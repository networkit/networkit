/*
 * RBMatrixReader.hpp
 *
 *  Created on: 16.10.2024
 *      Author: bernlu
 */

#ifndef NETWORKIT_IO_RB_MATRIX_READER_HPP_
#define NETWORKIT_IO_RB_MATRIX_READER_HPP_

#include <networkit/io/MatrixReader.hpp>
#include <networkit/io/RBGraphReader.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Reader for the Rutherford Boeing (RB) matrix file format as described in
 * http://sparse-files.engr.tamu.edu/files/DOC/rb.pdf.
 *
 * @note currently the reader only supports compressed column format for real, integer, or pattern
 * data types.
 */
class RBMatrixReader final : public MatrixReader {
public:
    RBMatrixReader() = default;

    CSRMatrix read(std::string_view path) override;

    /** Reads the matrix in @a in. */
    CSRMatrix read(std::istream &in);

private:
    friend RBGraphReader;
    // data format is CSC
    std::vector<index> pointers;
    std::vector<index> rowindex;
    std::vector<double> values;

    count nMatrixRows;
    count nMatrixCols;
    count nMatrixVals;
    bool symmetric = false;
    bool patternOnly = false;
    void readToVectors(std::istream &in);
};

} // namespace NetworKit
#endif // NETWORKIT_IO_RB_MATRIX_READER_HPP_

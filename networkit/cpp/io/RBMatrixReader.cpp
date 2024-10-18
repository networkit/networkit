/*
 *  RBMatrixReader.cpp
 *
 *  Created on: 16.10.2024
 *      Author: bernlu
 */

#include <fstream>
#include <stdexcept>
#include <string>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/io/RBMatrixReader.hpp>

namespace NetworKit {

namespace {

std::string tolower(std::string_view str) {
    std::string out;
    std::transform(str.begin(), str.end(), std::back_inserter(out), ::tolower);
    return out;
}
} // namespace

CSRMatrix RBMatrixReader::read(std::string_view path) {
    std::ifstream in(path.data());
    if (!in.is_open()) {
        throw std::runtime_error("could not open: " + std::string(path));
    }
    return read(in);
}

void RBMatrixReader::readToVectors(std::istream &in) {
    enum { POINTERS, ROWINDICES, VALUES } state = POINTERS;

    count totalRows;
    count nPointerRows;
    count nIndexRows;
    count nValueRows;
    std::string fmt;

    std::string line;
    std::istringstream in_line;

    // read header

    // first line
    getline(in, line);
    // pass - contains only text metadata

    // second line
    getline(in, line);
    in_line = std::istringstream(line);
    in_line >> totalRows >> nPointerRows >> nIndexRows >> nValueRows;

    // third line
    getline(in, line);
    in_line = std::istringstream(line);
    in_line >> fmt >> nMatrixCols >> nMatrixRows >> nMatrixVals;

    // make sure that the format is supported (first 3 characters are format information)
    if (tolower(fmt)[0] != 'r' && tolower(fmt)[0] != 'i' && tolower(fmt)[0] != 'p') {
        throw std::runtime_error(
            "Unsupported format: only real, integer, and pattern formats are supported.");
    }

    if (tolower(fmt)[0] == 'p')
        patternOnly = true;

    if (tolower(fmt)[1] == 's')
        symmetric = true;

    if (tolower(fmt)[2] != 'a') {
        throw std::runtime_error("Unsupported format: only compressed column format is supported.");
    }

    // fourth line
    getline(in, line);
    // pass - fortran pointer format is not relevant for this reader

    // prepare csc vectors
    pointers.reserve(nMatrixCols);
    rowindex.reserve(nMatrixVals);
    if (!patternOnly)
        values.reserve(nMatrixVals);
    else
        assert(nValueRows == 0); // there should be no value rows for p format.

    count row = 0;
    // start reading data
    while (getline(in, line)) {
        row++;
        std::istringstream in_line(line);

        if (state == POINTERS) {
            int p;
            while (in_line >> p) {
                pointers.push_back(
                    p - 1); // data in file is 1-indexed. shift to 0-index for networkit
            }

            if (row == nPointerRows) {
                state = ROWINDICES;
                row = 0;
            }
        }

        if (state == ROWINDICES) {
            int r;
            while (in_line >> r) {

                rowindex.push_back(
                    r - 1); // data in file is 1-indexed. shift to 0-index for networkit
            }

            if (row == nIndexRows) {
                state = VALUES;
                row = 0;
            }
        }

        if (state == VALUES) {
            double v;
            while (in_line >> v) {
                values.push_back(v);
            }

            if (row == nValueRows) {
                break;
            }
        }
    }
}

CSRMatrix RBMatrixReader::read(std::istream &in) {
    readToVectors(in);
    CSRMatrix mat(nMatrixRows, nMatrixCols, pointers, rowindex, values);

    // symmetric matrices are stored as triangular matrix.
    // Turn the data into a complete matrix again:
    if (symmetric)
        return mat + mat.transpose() - CSRMatrix::diagonalMatrix(mat.diagonal());
    // otherwise, transpose the matrix to convert from csc to csr.
    return mat.transpose();
}

} // namespace NetworKit

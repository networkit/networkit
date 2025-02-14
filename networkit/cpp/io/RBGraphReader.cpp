/*
 *  RBGraphReader.cpp
 *
 *  Created on: 16.10.2024
 *      Author: bernlu
 */

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/StringTools.hpp>
#include <networkit/io/RBGraphReader.hpp>
#include <networkit/io/RBMatrixReader.hpp>

namespace NetworKit {

Graph RBGraphReader::read(std::string_view path) {
    std::ifstream in(path.data());
    if (!in.is_open()) {
        throw std::runtime_error("could not open: " + std::string(path));
    }

    RBMatrixReader reader;

    reader.readToVectors(in);

    if (reader.nMatrixCols != reader.nMatrixRows)
        throw std::runtime_error(
            "File does not contain a square matrix - cannot parse this file into a graph!");

    Graph graph(reader.nMatrixCols, !reader.patternOnly, !reader.symmetric);

    // iterate values in csc matrix and add edges to graph
    for (index col = 0; col < reader.nMatrixCols; ++col) {
        for (index idx = reader.pointers[col]; idx <= reader.pointers[col + 1] - 1; ++idx) {
            if (reader.patternOnly)
                graph.addEdge(reader.rowindex[idx], col);
            else
                graph.addEdge(reader.rowindex[idx], col, reader.values[idx]);
        }
    }

    return graph;
}

} /* namespace NetworKit */

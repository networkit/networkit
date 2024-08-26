/*
 * DotWriter.cpp
 *
 *  Created on: Jun 5, 2013
 *      Author: forigem
 */

#include <fstream>

#include <networkit/io/DotGraphWriter.hpp>

namespace NetworKit {

void DotGraphWriter::write(const Graph &G, std::string_view path) {
    std::ofstream file{path.data()};

    file << "graph {\n";
    G.forEdges([&](node u, node v) { file << u << " -- " << v << ";\n"; });
    file << "}\n";
}

} /* namespace NetworKit */

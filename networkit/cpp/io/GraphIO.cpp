/*
 * GraphIO.cpp
 *
 *  Created on: 09.01.2013
 *      Author: Christian Staudt
 */

// networkit-format

#include <fstream>
#include <networkit/io/GraphIO.hpp>

namespace NetworKit {

void GraphIO::writeEdgeList(const Graph &G, const std::string &path) {

    std::ofstream file;
    file.open(path.c_str());

    G.forEdges([&](const node v, const node w) { file << v << " " << w << std::endl; });

    file.close();
    INFO("wrote graph to file: ", path);
}

void GraphIO::writeAdjacencyList(const Graph &G, const std::string &path) {
    std::ofstream file;
    file.open(path.c_str());

    G.forNodes([&](const node v) {
        file << v;
        G.forNeighborsOf(v, [&](const node x) { file << " " << x; });
        file << std::endl;
    });
    INFO("wrote graph to file: ", path);
}

} /* namespace NetworKit */

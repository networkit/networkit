/*
 * SNAPGraphWriter.cpp
 *
 *  Created on: 24.09.2013
 *      Author: cls
 */

#include "SNAPGraphWriter.h"
#include "../auxiliary/Enforce.h"

namespace NetworKit {

void SNAPGraphWriter::write(const Graph& G, const std::string& path) {
    std::ofstream file(path);
    Aux::enforceOpened(file);

    // write "problem line" - n, m, directed/undirected, weighted/weight type
    file << "p " << G.numberOfNodes() << " " << G.numberOfEdges() << " u u 0\n"; // FIXME: makeshift

    G.forEdges([&](node u, node v){
        file << u << " " << v << '\n';
    });

    file.close();
}

} /* namespace NetworKit */

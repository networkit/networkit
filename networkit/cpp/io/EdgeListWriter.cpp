/*
 * EdgeListWriter.cpp
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#include "EdgeListWriter.h"
#include "../auxiliary/Log.h"

#include <sstream>

#include "../auxiliary/Enforce.h"

namespace NetworKit {

EdgeListWriter::EdgeListWriter(char separator, node firstNode) : separator(separator), firstNode(firstNode) {}

void EdgeListWriter::write(const Graph& G, std::string path) {
    std::ofstream file(path);
    Aux::enforceOpened(file);

    if (G.isWeighted()) {
        G.forEdges([&](node u, node v, double weight){
            file << (u + firstNode) << separator << (v + firstNode) << separator << weight << std::endl;
        });
    } else {
        G.forEdges([&](node u, node v){
        	file << (u + firstNode) << separator << (v + firstNode) << std::endl;
        });
    }


    file.close();

}

} /* namespace NetworKit */

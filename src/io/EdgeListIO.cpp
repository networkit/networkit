/*
 * EdgeListIO.cpp
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#include "EdgeListIO.h"

namespace NetworKit {


EdgeListIO::EdgeListIO() {

}

EdgeListIO::EdgeListIO(char separator, node firstNode) : separator(separator), firstNode(firstNode) {

}

EdgeListIO::~EdgeListIO() {

}



Graph EdgeListIO::read(std::string path) {

    std::ifstream file;
    std::string line; // the current line

    // read file once to get to the last line and figure out the number of nodes
    // unfortunately there is an empty line at the ending of the file, so we need to get the line before that

   file.open(path);

   std::string previousLine;
   node maxNode = 0;

   std::string commentPrefix = "#";

   // first find out the maximum node id
   DEBUG("first pass");
   while (file.good()) {

        std::getline(file, line);

        if (line.compare(0, commentPrefix.length(), commentPrefix) == 0) {
            // TRACE("ignoring comment: " << line);
        } else if (line.length() == 0) {
            // TRACE("ignoring empty line");
        } else {
            // TRACE("edge line: " << line);
            std::vector<std::string> split = Aux::StringTools::split(line, separator);
       
            if (split.size() == 2) {

                node u = std::stoi(split[0]);
                if (u > maxNode) {
                    maxNode = u;
                }
                node v = std::stoi(split[1]);
                if (v > maxNode) {
                    maxNode = v;
                }
            } else {
                throw std::runtime_error("malformed line in edge list file");
            }
        }

    }

    maxNode = maxNode - firstNode + 1;
    TRACE("max. node id found: " << maxNode);

    Graph G(maxNode);

    file.close();

    DEBUG("second pass");
    file.open(path);

    // split the line into start and end node. since the edges are sorted, the start node has the highest id of all nodes
    while(std::getline(file,line)){

        if (line.compare(0, commentPrefix.length(), commentPrefix) == 0) {
            // TRACE("ignoring comment: " << line);
        } else {
            // TRACE("edge line: " << line);
            std::vector<std::string> split = Aux::StringTools::split(line, separator);
            std::string splitZero = split[0];
            if (split.size() == 2) {
                node u = std::stoi(split[0]) - firstNode;
                node v = std::stoi(split[1]) - firstNode;
                if (!G.hasEdge(u,v) && !G.hasEdge(v,u)) {
                    G.addEdge(u, v);
                }
            } else {
                throw std::runtime_error("malformed line in edge list file");
            }
        }
    }

    file.close();

    return G;
}



void EdgeListIO::write(const Graph& G, std::string path) {
    std::ofstream file;
    file.open(path);

    assert (file.good());

    G.forEdges([&](node u, node v){
    	file << (u + firstNode) << " " << (v + firstNode) << std::endl;
    });

    file.close();

}

} /* namespace NetworKit */

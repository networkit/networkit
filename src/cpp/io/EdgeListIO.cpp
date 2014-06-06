/*
 * EdgeListIO.cpp
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#include "EdgeListIO.h"
#include "../auxiliary/Log.h"

#include <sstream>

#include "../auxiliary/Enforce.h"

namespace NetworKit {

EdgeListIO::EdgeListIO(char separator, node firstNode) : separator(separator), firstNode(firstNode) {}

Graph EdgeListIO::read(const std::string& path) {

    std::ifstream file(path);
    Aux::enforceOpened(file);
    std::string line; // the current line

    // read file once to get to the last line and figure out the number of nodes
    // unfortunately there is an empty line at the ending of the file, so we need to get the line before that

   std::string previousLine;
   node maxNode = 0;

   std::string commentPrefix = "#";

   DEBUG("separator: " , separator);
   DEBUG("first node: " , firstNode);
   // first find out the maximum node id
   DEBUG("first pass");
   count i = 0;
   while (file.good()) {
        ++i;
        std::getline(file, line);
        // TRACE("read line: " , line);

        if (line.compare(0, commentPrefix.length(), commentPrefix) == 0) {
            // TRACE("ignoring comment: " , line);
        } else if (line.length() == 0) {
            // TRACE("ignoring empty line");
        } else {
            std::vector<std::string> split = Aux::StringTools::split(line, separator);
       
            if (split.size() == 2) {
                TRACE("split into : " , split[0] , " and " , split[1]);
                node u = std::stoi(split[0]);
                if (u > maxNode) {
                    maxNode = u;
                }
                node v = std::stoi(split[1]);
                if (v > maxNode) {
                    maxNode = v;
                }
            } else {
                std::stringstream message;
                message << "malformed line ";
                message << i << ": ";
                message << line;
                throw std::runtime_error(message.str());
            }
        }

    }

    file.close();

    maxNode = maxNode - firstNode + 1;
    DEBUG("max. node id found: " , maxNode);

    Graph G(maxNode);



    DEBUG("second pass");
    file.open(path);

    // split the line into start and end node. since the edges are sorted, the start node has the highest id of all nodes
    i = 0; // count lines
    while(std::getline(file,line)){
        ++i;
        if (line.compare(0, commentPrefix.length(), commentPrefix) == 0) {
            // TRACE("ignoring comment: " , line);
        } else {
            // TRACE("edge line: " , line);
            std::vector<std::string> split = Aux::StringTools::split(line, separator);
            std::string splitZero = split[0];
            if (split.size() == 2) {
                node u = std::stoi(split[0]) - firstNode;
                node v = std::stoi(split[1]) - firstNode;
                if (!G.hasEdge(u,v) && !G.hasEdge(v,u)) {
                    G.addEdge(u, v);
                }
            } else {
                std::stringstream message;
                message << "malformed line ";
                message << i << ": ";
                message << line;
                throw std::runtime_error(message.str());
            }
        }
    }

    file.close();

    return G;
}



void EdgeListIO::write(const Graph& G, std::string path) {
    std::ofstream file(path);
    Aux::enforceOpened(file);

    G.forEdges([&](node u, node v){
    	file << (u + firstNode) << separator << (v + firstNode) << std::endl;
    });

    file.close();

}

} /* namespace NetworKit */

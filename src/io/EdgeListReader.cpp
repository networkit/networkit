/*
 * EdgeListReader.cpp
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#include "EdgeListReader.h"

namespace NetworKit {

EdgeListReader::EdgeListReader(node firstNode) : firstNode(firstNode) {
	// TODO Auto-generated constructor stub

}

EdgeListReader::~EdgeListReader() {
	// TODO Auto-generated destructor stub
}



Graph EdgeListReader::read(std::string path) {

    std::ifstream file;
    std::string line; // the current line

    // read file once to get to the last line and figure out the number of nodes
    // unfortunately there is an empty line at the ending of the file, so we need to get the line before that

   file.open(path);

   std::string previousLine;
   node maxNode = 0;

   while (file.good()) {
	   std::getline(file, line);
	   std::vector<std::string> split = Aux::StringTools::split(line, '\t');
	   std::string prefix = "#";
	   if (split.size() == 2 && (split[0].compare(0, prefix.length(), prefix) != 0)) {
		   node u = std::stoi(split[0]);
		   if (u > maxNode) {
			   maxNode = u;
		   }
		   node v = std::stoi(split[1]);
		   if (v > maxNode) {
			   maxNode = v;
		   }
	   }
    }

    maxNode = maxNode - firstNode + 1;
    Graph G(maxNode);

    file.close();

    file.open(path);

    // split the line into start and end node. since the edges are sorted, the start node has the highest id of all nodes
    while(std::getline(file,line)){
    	std::vector<std::string> split = Aux::StringTools::split(line, '\t');
    	std::string prefix = "#";
    	std::string splitZero = split[0];
    	bool notHash = (splitZero.compare(0, prefix.length(), prefix) != 0);
		if (split.size() == 2 && notHash) {
			node u = std::stoi(split[0]) - firstNode;
			node v = std::stoi(split[1]) - firstNode;
			if (!G.hasEdge(u,v) && !G.hasEdge(v,u)) {
				G.addEdge(u, v);
			}
		}
    }

    file.close();



    return G;
}

} /* namespace NetworKit */

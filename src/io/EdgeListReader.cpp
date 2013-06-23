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
    while (file.good()) {
    	previousLine = line;
    	std::getline(file, line);
    }

    std::vector<std::string> split = Aux::StringTools::split(previousLine, '\t');
    node u = std::stoi(split[0]) - firstNode + 1;
    Graph G(u);




    file.close();

    file.open(path);
    // Graph G(0);

    // split the line into start and end node. since the edges are sorted, the start node has the highest id of all nodes
    while(std::getline(file,line)){
    	std::vector<std::string> split = Aux::StringTools::split(line, '\t');
		if (split.size() == 2 && split[0] != "#") {
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

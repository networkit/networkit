/*
 * EdgeListReader.cpp
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#include "EdgeListReader.h"

namespace NetworKit {

EdgeListReader::EdgeListReader() {
	// TODO Auto-generated constructor stub

}

EdgeListReader::~EdgeListReader() {
	// TODO Auto-generated destructor stub
}

Graph EdgeListReader::read(std::string path) {
    std::ifstream file;
    std::string line; // the current line

    // read file once to get to the last line and figure out the number of nodes
    file.open(path);

    while (file.good()) {
    	std::getline(file, line);
    }
    // TODO: unfortunately there is an empty line at the ending of the file, so we need to get the line before that
    DEBUG("the last line is: " << line);

    // TODO: split the line into start and end node. since the edges are sorted, the start node has the highest id of all nodes

    file.close();

    file.open(path);
    // TODO: read the file line by line and insert edges
    // TODO: remember that nodes in this format is 1-based while NetworKit nodes are 0-based. the base index could also be a parameter of the class
    // to make it more flexible
    file.close();


    Graph G(0);
    return G;
}

} /* namespace NetworKit */

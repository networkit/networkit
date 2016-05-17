/*
 * KONECTGraphReader.cpp
 * 
 * Reader for the KONECT graph format, 
 * based on the EdgeListReader.
 * 
 * The KONECT format is described in detail in 
 * http://konect.uni-koblenz.de/downloads/konect-handbook.pdf
 */

#include "KONECTGraphReader.h"
#include "../auxiliary/Log.h"

#include <sstream>

#include "../auxiliary/Enforce.h"

#include <algorithm>

namespace NetworKit {

    KONECTGraphReader::KONECTGraphReader(char separator, bool ignoreLoops) :
    separator(separator), commentPrefix("%"), firstNode(1), ignoreLoops(ignoreLoops) {
    }

    Graph KONECTGraphReader::read(const std::string& path) {
        return readContinuous(path);
    }
    Graph KONECTGraphReader::readContinuous(const std::string& path) {

        std::ifstream file(path);
        Aux::enforceOpened(file);
        std::string line; // the current line


        DEBUG("separator: ", this->separator);
        DEBUG("first node: ", this->firstNode);

        // first find out the maximum node id
        DEBUG("first pass");
 	// first run through the file determines if the graph is directed and/or weighted and checks the consistency of the file
	// IF NEEDED: try to improve performance by storing edges in a vector or map during first pass.
	node maxNode = 0;
        count i = 0;

        // the number of vertices / edges as specified in the input file
        int n = -1, m = -1;
	// directed or weighted graph?
        bool directed = true, weighted = false;
        // the minimum number of values per line
        unsigned int minValuesPerLine = 2;
	// attempt to detect a tab separator character 
	bool detectSeparatorAttempt = true;


        while (file.good()) {
            ++i;
            std::getline(file, line);
            // TRACE("read line: " , line);
            if (line.compare(0, this->commentPrefix.length(), this->commentPrefix) == 0) {
                if (i == 1) {
		    // first comment line determines if graph is directed/undirected and weighted/unweighted
                    std::vector<std::string> split = Aux::StringTools::split(line, ' ');
                    if (split.size() < 3) {
                        std::stringstream message;
                        message << "malformed line - first line must contain graph format and weight information, in line ";
                        message << i << ": " << line;
                        throw std::runtime_error(message.str());
                    }
                    if (split[1] == "sym" || split[1] == "bip") {
                        directed = false;
                    } else if (split[1] == "asym") {
                        directed = true;
                    } else {
                        std::stringstream message;
                        message << "malformed line - first line must give the graph format (sym/asym/bip), in line ";
                        message << i << ": " << line;
                        throw std::runtime_error(message.str());
                    }

                    if (split[2] == "unweighted" || split[2] == "positive") {
                        weighted = false;
			// NOTE: positive means, that there can be multiple edges. currently, these will be ignored. 
                        // graph must only contain source and target ids
                        minValuesPerLine = 2;
                    } else if (split[2] == "posweighted" || split[2] == "signed" || split[2] == "weighted") {
                        weighted = true;
                        // graph must contain source and target ids and weight!
                        minValuesPerLine = 3;
                    } else {
			std::stringstream message;
			message << "graph types \"multiweighted\" and \"dynamic\" are not supported yet, found in line ";
			message << i << ": " << line;
			throw std::runtime_error(message.str());
		    }
                } else if (i == 2) {
		    // the second, optional comment line contains number of vertices / edges
                    std::vector<std::string> split = Aux::StringTools::split(line, ' ');
                    m = std::stoul(split[1]);
                    n = std::stoul(split[2]);
//		    TRACE("m is: ",m);
//		    TRACE("n is: ",n);
                }
            } else if (line.length() == 0) {
                // TRACE("ignoring empty line");
            } else {
                std::vector<std::string> split = Aux::StringTools::split(line, this->separator);
                split.erase(std::remove_if(split.begin(),split.end(),[](const std::string& s){return s.empty();}),split.end());

		if(detectSeparatorAttempt) {
                    // one attempt is made to detect if, instead of a space, 
                    // a tab separator is used in the input file
		    detectSeparatorAttempt = false;
		    if(separator == ' ') {
                	std::vector<std::string> tabSplit = Aux::StringTools::split(line, '\t');
                        if(tabSplit.size() >= 2 && tabSplit.size() > split.size()) {
                            DEBUG("detected tab separator.");
                            split = tabSplit;
                            this->separator = '\t';
                        }
		    }
                }


                // correct number of values?
                if (split.size() >= minValuesPerLine && split.size() <= 4) {
                    TRACE("split into : ", split[0], " and ", split[1]);
                    node u = std::stoul(split[0]);
                    node v = std::stoul(split[1]);
                    if (!this->ignoreLoops || u != v) {
                        if (u > maxNode) {
                            maxNode = u;
                        }
                        if (v > maxNode) {
                            maxNode = v;
                        }
                    }
                } else {
                    std::stringstream message;
                    message << "malformed line (expecting ";
                    message << minValuesPerLine << "-4 values, ";
                    message << split.size() << " given) ";
                    message << i << ": ";
                    for (const auto& s : split) {
                        message << s <<", ";
                    }
                    throw std::runtime_error(message.str());
                }
            }
        }
        file.close();
        maxNode = maxNode - this->firstNode + 1;
        DEBUG("max. node id found: ", maxNode);

        Graph G(maxNode, weighted, directed);

        DEBUG("second pass");
	// second pass adds the edges to the graph.
        file.open(path);
        i = 0; // count lines
        while (std::getline(file, line)) {
            ++i;
            if (line.compare(0, this->commentPrefix.length(), this->commentPrefix) == 0) {
                // TRACE("ignoring comment: " , line);
                // comment lines already processed in first pass
            } else {
                // TRACE("edge line: " , line);
                std::vector<std::string> split = Aux::StringTools::split(line, this->separator);
                split.erase(std::remove_if(split.begin(),split.end(),[](const std::string& s){return s.empty();}),split.end());
                // correct number of values?
                if (split.size() >= minValuesPerLine && split.size() <= 4) {
                    node u = std::stoul(split[0]) - this->firstNode;
                    node v = std::stoul(split[1]) - this->firstNode;
                    if (weighted) {
                        count weightIdx = (split[2].size() > 0) ? 2 : 3;
                        edgeweight weight = std::stod(split[weightIdx]);

                        if (!this->ignoreLoops || u != v) {
                            if (directed) {
                                if (!G.hasEdge(u, v)) {
                                    G.addEdge(u, v, weight);
                                }
                            } else {
                                if (!G.hasEdge(u, v) && !G.hasEdge(v, u)) {
                                    G.addEdge(u, v, weight);
                                }
                            }
                        }
                    } else {

                        if (!this->ignoreLoops || u != v) {
                            if (directed) {
                                if (!G.hasEdge(u, v)) {
                                    G.addEdge(u, v);
                                }
                            } else {
                                if (!G.hasEdge(u, v) && !G.hasEdge(v, u)) {
                                    G.addEdge(u, v);
                                }
                            }
                        }
                    }
                } else {
                    std::stringstream message;
                    message << "malformed line (expecting ";
                    message << minValuesPerLine << "-4 values, ";
                    message << split.size() << " given) ";
                    message << i << ": " << line;
                    throw std::runtime_error(message.str());
                }
            }
        }
        file.close();
        if (m != -1 && G.numberOfEdges() != (unsigned int) m) {
            ERROR("KONECT file is corrupted: actual number of added edges doesn't match the specifed number of edges");
        }
        if (n != -1 && G.numberOfNodes() != (unsigned int) n) {
            ERROR("KONECT file is corrupted: actual number of added vertices doesn't match the specifed number of vertices");
        }

        G.shrinkToFit();
        return G;
    }


} /* namespace NetworKit */

/*
 * METISParser.cpp
 *
 *  Created on: 27.11.2012
 *      Author: Christian Staudt
 */

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/io/METISParser.hpp>

#include <stdexcept>

namespace NetworKit {

METISParser::METISParser(std::string_view path) : graphFile(path.data()) {
    if (!(this->graphFile)) {
        ERROR("invalid graph file: ", path);
        throw std::runtime_error("invalid graph file");
    }
}

std::tuple<count, count, index, count> METISParser::getHeader() {
    // handle header line
    count n;        // number of nodes
    count m;        // number of edges
    index fmt = 0;  // weighted or unweighted graph
    count ncon = 0; // number of node weights

    std::string line = "";
    Aux::enforceOpened(this->graphFile);

    if (std::getline(this->graphFile, line)) {
        // ignore comment lines
        while (line[0] == '%') {
            std::getline(this->graphFile, line);
        }

        std::vector<count> tokens = METISParserDetails::parseLine<count>(line);
        n = tokens[0];
        m = tokens[1];
        if (tokens.size() == 2) {
            return std::tuple<count, count, index, count>(n, m, fmt, ncon);
        }
        if (tokens.size() >= 3) {
            fmt = tokens[2];
            if (fmt >= 2) {
                WARN("nodes are weighted; node weights will be ignored");
            }
            if (tokens.size() == 4) {
                ncon = tokens[3];
            } else {
                ncon = 1;
            }
        }
        return std::tuple<count, count, index, count>(n, m, fmt, ncon);
    } else {
        ERROR("getline not successful");
        throw std::runtime_error("getting METIS file header failed");
        return std::tuple<count, count, index, count>(0, 0, 0, 0);
    }
}

bool METISParser::hasNext() {
    // if graph file has lines left, return true
    return this->graphFile.good();
}

std::string METISParser::nextLine() {
    std::string line;
    do {
        std::getline(this->graphFile, line);
        // check for comment line starting with '%'
    } while (line[0] == '%');
    return line;
}

} /* namespace NetworKit */

/*
 * DGSStreamParser.hpp
 *
 *  Created on: 23.12.2013
 *      Author: cls
 */

#ifndef NETWORKIT_DYNAMICS_DGS_STREAM_PARSER_HPP_
#define NETWORKIT_DYNAMICS_DGS_STREAM_PARSER_HPP_

#include <fstream>
#include <map>
#include <string>
#include <vector>

#include <networkit/dynamics/GraphEvent.hpp>

namespace NetworKit {

/**
 * @ingroup dynamics
 */
class DGSStreamParser final {

public:

    DGSStreamParser(std::string path, bool mapped = true, node baseIndex = 0);

    std::vector<GraphEvent> getStream();

private:

    std::ifstream dgsFile;
    bool mapped;
    std::map<std::string, node> key2id;
    node baseIndex;
    node nextNode;
};

} /* namespace NetworKit */

#endif // NETWORKIT_DYNAMICS_DGS_STREAM_PARSER_HPP_

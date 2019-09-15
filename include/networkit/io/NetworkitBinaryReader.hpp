/*
 * NetworkitBinaryReader.hpp
 *
 *      Author: Charmaine Ndolo <charmaine.ndolo@b-tu.de>
 */

#ifndef NETWORKIT_BINARY_READER_H
#define NETWORKIT_BINARY_READER_H

#include <networkit/io/GraphReader.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup io
 *
 * Reads a graph written in the custom Networkit format documented in cpp/io/NetworkitGraph.md
 */

class NetworkitBinaryReader : public GraphReader {

public:
	NetworkitBinaryReader() {};

	Graph read(const std::string& path) override; 

private:
	count nodes;
	count chunks;
	bool directed;
	bool weighted;
};
} /* namespace */

#endif

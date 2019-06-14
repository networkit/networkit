/*
 * NetworkitBinaryWriter.h
 *
 *@author Charmaine Ndolo <charmaine.ndolo@b-tu.de>
 */

#ifndef NETWORKIT_BINARY_WRITER_H
#define NETWORKIT_BINARY_WRITER_H

#include <networkit/graph/Graph.hpp>
#include <networkit/io/GraphWriter.hpp>
#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/io/NetworkitBinaryGraph.hpp>
#include <tlx/math/clz.hpp>

#include <fstream>
#include <cstring>

namespace NetworKit {

/**
 * @ingroup io
 *
 * Writes a graph in the custom Networkit format documented in cpp/io/NetworkitGraph.md
 */
	
class NetworkitBinaryWriter : public GraphWriter {
	
public:
	NetworkitBinaryWriter(uint64_t chunks = 32);
	
	void write(const Graph& G, const std::string& path);
		
private:
	static size_t encode(uint64_t value, uint8_t* buffer);

	static uint64_t encodeZigzag(int64_t value);

	count nodes;
	count chunks;	
};
} /* namespace */
#endif

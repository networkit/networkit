/*
 * NetworkitBinaryWriter.hpp
 *
 * @author Charmaine Ndolo <charmaine.ndolo@b-tu.de>
 */

#ifndef NETWORKIT_BINARY_WRITER_H
#define NETWORKIT_BINARY_WRITER_H

#include <networkit/graph/Graph.hpp>
#include <networkit/io/GraphWriter.hpp>

namespace NetworKit {

/**
 * @ingroup io
 *
 * Writes a graph in the custom Networkit format documented in cpp/io/NetworkitGraph.md
 */
enum class NetworkitBinaryWeights {
	none,
	unsignedFormat,
	signedFormat,
	doubleFormat,
	floatFormat,
	autoDetect
};

class NetworkitBinaryWriter final : public GraphWriter {

public:
	NetworkitBinaryWriter(uint64_t chunks = 32, NetworkitBinaryWeights weightsType = NetworkitBinaryWeights::autoDetect);

	void write(const Graph &G, const std::string &path) override;

private:
	static size_t encode(uint64_t value, uint8_t* buffer);

	static uint64_t encodeZigzag(int64_t value);

	count chunks;
	NetworkitBinaryWeights weightsType;
};
} /* namespace */
#endif

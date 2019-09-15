/*
 * NetworkitBinaryWriter.hpp
 *
 *      Author: Charmaine Ndolo <charmaine.ndolo@b-tu.de>
 */

#ifndef NETWORKIT_BINARY_WRITER_H
#define NETWORKIT_BINARY_WRITER_H

#include <networkit/graph/Graph.hpp>
#include <networkit/io/GraphWriter.hpp>

namespace NetworKit {

enum class NetworkitBinaryWeights {
	none,
	unsignedFormat,
	signedFormat,
	doubleFormat,
	floatFormat,
	autoDetect
};

/**
 * @ingroup io
 *
 * Writes a graph in the custom Networkit format documented in cpp/io/NetworkitGraph.md
 */
class NetworkitBinaryWriter final : public GraphWriter {

public:
	NetworkitBinaryWriter(uint64_t chunks = 32, NetworkitBinaryWeights weightsType = NetworkitBinaryWeights::autoDetect);

	void write(const Graph &G, const std::string &path) override;

private:
	count chunks;
	NetworkitBinaryWeights weightsType;
};

} // namespace NetworKit

#endif

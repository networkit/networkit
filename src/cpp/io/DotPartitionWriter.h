/*
 * DotPartitionWriter.h
 */

#ifndef DOTPARTITIONWRITER_H_
#define DOTPARTITIONWRITER_H_

#include <map>

#include "../structures/Partition.h"
#include "../graph/Graph.h"

namespace NetworKit {

class DotPartitionWriter {
public:
	virtual ~DotPartitionWriter() = default;
	
	virtual void write(Graph& graph, Partition& zeta, const std::string& path) const;
	
private:
	std::size_t countConnectedNodes(const Graph& graph) const;

};

} /* namespace NetworKit */
#endif /* DOTPARTITIONGWRITER_H_ */

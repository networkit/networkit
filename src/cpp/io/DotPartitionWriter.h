/*
 * DotPartitionWriter.h
 */

#ifndef DOTPARTITIONWRITER_H_
#define DOTPARTITIONWRITER_H_

#include <fstream>

#include "../structures/Partition.h"
#include "../graph/Graph.h"

namespace NetworKit {

class DotPartitionWriter {

public:
    DotPartitionWriter();

    virtual ~DotPartitionWriter();

    virtual void write(Graph& graph, Partition& zeta, std::string path) const;

    virtual std::map<index, double> createHueMap(Graph &graph, Partition& zeta) const;
};

} /* namespace NetworKit */
#endif /* DOTPARTITIONGWRITER_H_ */

/*
 * DotClusteringWriter.h
 */

#ifndef DOTCLUSTERINGWRITER_H_
#define DOTCLUSTERINGWRITER_H_

#include <fstream>

#include "../structures/Partition.h"
#include "../graph/Graph.h"

namespace NetworKit {

class DotClusteringWriter {

public:
    DotClusteringWriter();

    virtual ~DotClusteringWriter();

    virtual void write(Graph& graph, Partition& zeta, std::string path) const;

    virtual std::map<index, double> createHueMap(Graph &graph, Partition& zeta) const;
};

} /* namespace NetworKit */
#endif /* CLUSTERINGWRITER_H_ */

/*
 * DotClusteringWriter.h
 */

#ifndef DOTCLUSTERINGWRITER_H_
#define DOTCLUSTERINGWRITER_H_

#include <fstream>

#include "../clustering/Clustering.h"

namespace NetworKit {

class DotClusteringWriter {

public:
    DotClusteringWriter();

    virtual ~DotClusteringWriter();

    virtual void write(Graph& graph, Clustering& zeta, std::string path) const;

    virtual std::map<cluster, double> createHueMap(Graph &graph, Clustering& zeta) const;
};

} /* namespace NetworKit */
#endif /* CLUSTERINGWRITER_H_ */

/*
 * DotClusteringWriter.h
 */

#ifndef DOTWRITER_H_
#define DOTWRITER_H_

#include <fstream>

#include "../clustering/Clustering.h"

namespace NetworKit {

class DotWriter {

public:
    DotWriter();

    virtual ~DotWriter();

    virtual void write(Graph& graph, std::string path) const;

};

} /* namespace NetworKit */
#endif /* DOTWRITER_H_ */

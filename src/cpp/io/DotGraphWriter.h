/*
 * DotGraphWriter.h
 */

#ifndef DOTGRAPHWRITER_H
#define DOTGRAPHWRITER_H

#include <fstream>
#include "../graph/Graph.h"

namespace NetworKit {

class DotGraphWriter {
public:

    virtual ~DotGraphWriter() = default;

    virtual void write(Graph& graph, std::string path) const;

};

} /* namespace NetworKit */
#endif /* DOTGRAPHWRITER_H */

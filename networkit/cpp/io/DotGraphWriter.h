/*
 * DotGraphWriter.h
 */

#ifndef DOTGRAPHWRITER_H
#define DOTGRAPHWRITER_H

#include <fstream>
#include "../graph/Graph.h"

namespace NetworKit {

/**

 * @ingroup io
 *
 * This class turns a graph into a very basic GraphViz file as documented in the official manual [1].
 * If a more thorough support is desired, please contact the developers over networkit@ira.uni-karlsruhe.de.
 *
 * [1] http://www.graphviz.org/Documentation/dotguide.pdf
 */ 
class DotGraphWriter {
public:
	/**
	 * Write a graph as a GraphViz/file.
	 * 
	 * @param[in]	graph	The graph object
	 * @param[in]	path	The file path to be written to
	 */
	virtual void write(Graph& graph, std::string path) const;

};

} /* namespace NetworKit */
#endif /* DOTGRAPHWRITER_H */

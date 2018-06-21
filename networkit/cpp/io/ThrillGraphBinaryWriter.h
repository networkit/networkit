/*
 * ThrillGraphBinaryWriter.h
 *
 * @author Michael Hamann <michael.hamann@kit.edu>
 */

#ifndef THRILLGRAPHBINARYWRITER_H_
#define THRILLGRAPHBINARYWRITER_H_

#include "GraphWriter.h"

namespace NetworKit {

class ThrillGraphBinaryWriter : public GraphWriter {
public:
	/**
	 * Write the given graph into a binary file at the given path.
	 * 
	 * @param[in] G The graph to write.
	 * @param[in] path The path where to write the graph.
	 */
	virtual void write(const Graph& G, const std::string& path);
};

} /* namespace NetworKit */


#endif /* THRILLGRAPHBINARYWRITER_H_ */

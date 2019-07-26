/*
 * ThrillGraphBinaryWriter.hpp
 *
 * @author Michael Hamann <michael.hamann@kit.edu>
 */

#ifndef THRILLGRAPHBINARYWRITER_H_
#define THRILLGRAPHBINARYWRITER_H_

#include <networkit/io/GraphWriter.hpp>

namespace NetworKit {

class ThrillGraphBinaryWriter final : public GraphWriter {
public:
	/**
	 * Write the given graph into a binary file at the given path.
	 *
	 * @param[in] G The graph to write.
	 * @param[in] path The path where to write the graph.
	 */
	void write(const Graph &G, const std::string &path) override;
};

} /* namespace NetworKit */


#endif /* THRILLGRAPHBINARYWRITER_H_ */

/*
 * CoreDecomposition.h
 *
 *  Created on: Oct 28, 2013
 *      Author: Henning
 */

#ifndef COREDECOMPOSITION_H_
#define COREDECOMPOSITION_H_

#include <vector>
#include <fstream>
#include <string>
#include "../graph/Graph.h"
#include "../auxiliary/ShellList.h"

namespace NetworKit {

/**
 * Computes k-core decomposition of a graph.
 */
class CoreDecomposition {

private:
	std::vector<count> coreness;

public:
	CoreDecomposition();
	virtual ~CoreDecomposition();

	/**
	 * @return k-core decomposition of graph @a G.
	 */
	std::vector<count> run(const Graph& G);

	void write(std::string filename);
};

} /* namespace NetworKit */
#endif /* COREDECOMPOSITION_H_ */

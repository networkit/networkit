/*
 * DynBFS.h
 *
 *  Created on: 17.07.2014
 *      Author: cls, ebergamini
 */

#ifndef DYNBFS_H_
#define DYNBFS_H_

#include "DynSSSP.h"
#include "BFS.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Dynamic breadth-first search.
 */
class DynBFS : public DynSSSP, public BFS {

public:

	/**
	 * Creates the object for @a G and source @a s.
	 *
	 * @param G The graph.
	 * @param s The source node.
	 */
	DynBFS(const Graph& G, node s);

	void init();

	/** Updates the distances after an event.*/
	void update(const std::vector<GraphEvent>& batch) override;

protected:
	std::vector<Color> color;
	count maxDistance;

};

} /* namespace NetworKit */

#endif /* DYNSSSP_H_ */

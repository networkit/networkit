/*
 * DynBFS.h
 *
 *  Created on: 17.07.2014
 *      Author: cls, ebergamini
 */

#ifndef DYNBFS_H_
#define DYNBFS_H_

#include "DynSSSP.h"

namespace NetworKit {

enum Color {WHITE, BLACK, GRAY};

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

	/** Default destructor */
	virtual ~DynBFS() = default;


    void update(const std::vector<GraphEvent>& batch) override;

protected:

	std::vector<Color> color;

};

} /* namespace NetworKit */

#endif /* DYNSSSP_H_ */

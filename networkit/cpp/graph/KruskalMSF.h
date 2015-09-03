/*
 * KruskalMSF.h
 *
 *  Created on: 03.09.2015
 *      Author: Henning
 */

#ifndef KRUSKALMSF_H_
#define KRUSKALMSF_H_

#include "../Globals.h"
#include "Graph.h"

namespace NetworKit {

class KruskalMSF {
protected:
	const Graph& G;
	Graph tree;

public:
	KruskalMSF(const Graph& G);
	virtual ~KruskalMSF() = default;

	void run();

	Graph getTree();
};

} /* namespace NetworKit */
#endif /* KRUSKALMSF_H_ */

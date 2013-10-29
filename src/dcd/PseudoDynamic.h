/*
 * PseudoDynamic.h
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#ifndef PSEUDODYNAMIC_H_
#define PSEUDODYNAMIC_H_

#include "../generators/DynamicGraphSource.h"

namespace NetworKit {

class PseudoDynamic: public NetworKit::DynamicGraphSource {
public:

	PseudoDynamic(Graph& Gstatic);

	virtual ~PseudoDynamic();

	virtual void initializeGraph();

	virtual void generate();


private:

	const Graph Gstatic; // cannot be a reference or a pointer because this object might survive longer than static graph

	std::vector<node> idmap; // maping from static graph id to dynamic graph id
	node current; // current node
};

} /* namespace NetworKit */
#endif /* PSEUDODYNAMIC_H_ */

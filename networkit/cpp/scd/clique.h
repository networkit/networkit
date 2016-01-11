#ifndef CLIQUE_H_
#define CLIQUE_H_

#include "../graph/Graph.h"

namespace NetworKit {

class Clique {

public:
	Clique(const Graph& G);

	std::vector<std::vector<node> > run(node& seed);

protected:
	const Graph& G;

	std::vector<std::vector<node> > tomita(std::vector<node>& pxvector, std::vector<node>& pxlookup, std::vector<std::vector<node> >& neighbors, uint32_t xbound, uint32_t xpbound, uint32_t pbound, std::vector<node>& r);

	node findPivot(std::vector<node>& pxvector, std::vector<node>& pxlookup, std::vector<std::vector<node> >& neighbors, uint32_t xbound, uint32_t xpbound, uint32_t pbound);

};

}

#endif /* CLIQUE_H_ */

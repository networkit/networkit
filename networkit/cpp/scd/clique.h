#ifndef CLIQUE_H_
#define CLIQUE_H_

#include "../graph/Graph.h"

#include <unordered_map>

namespace NetworKit {

class Clique {

public:
	Clique(const Graph& G);

	std::vector<std::vector<node> > run();

protected:
	const Graph& G;

	std::vector<std::vector<node> > tomita(std::vector<node>& pxvector, std::unordered_map<node, uint32_t>& pxlookup, std::vector<std::vector<node> >& neighbors, uint32_t xbound, uint32_t xpbound, uint32_t pbound, std::vector<node>& r);

	node findPivot(std::vector<node>& pxvector, std::unordered_map<node, uint32_t>& pxlookup, std::vector<std::vector<node> >& neighbors, uint32_t xbound, uint32_t xpbound, uint32_t pbound);

	std::vector<node> getDegeneracyOrdering();

};

}

#endif /* CLIQUE_H_ */

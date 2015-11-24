#ifndef CLIQUE_H_
#define CLIQUE_H_

#include "../graph/Graph.h"

#include <set>

namespace NetworKit {

class Clique {

public:
	Clique(const Graph& G, const unsigned int missingEdges);

	std::vector<std::set<node> > run(node& seed);

protected:
	const Graph& G;
	const unsigned int missingEdges;
};

}

#endif /* CLIQUE_H_ */

#ifndef CLIQUE_H_
#define CLIQUE_H_

#include "../graph/Graph.h"
#include "../base/Algorithm.h"

#include <unordered_map>

namespace NetworKit {

class MaximalCliques : public Algorithm {

public:
	MaximalCliques(const Graph& G);

	void run() override;

	const std::vector<std::vector<node>>& getCliques() const;

protected:
	const Graph& G;

	std::vector<std::vector<node>> result;
};

}

#endif /* CLIQUE_H_ */

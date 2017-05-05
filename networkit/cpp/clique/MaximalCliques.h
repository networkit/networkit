#ifndef CLIQUE_H_
#define CLIQUE_H_

#include "../graph/Graph.h"
#include "../base/Algorithm.h"

#include <unordered_map>

namespace NetworKit {

class MaximalCliques : public Algorithm {

public:
	MaximalCliques(const Graph& G);
	MaximalCliques(const Graph& G, std::function<void(const std::vector<node>&)> callback);

	void run() override;

	const std::vector<std::vector<node>>& getCliques() const;

protected:
	const Graph& G;

	std::vector<std::vector<node>> result;

	std::function<void(const std::vector<node>&)> callback;
};

}

#endif /* CLIQUE_H_ */

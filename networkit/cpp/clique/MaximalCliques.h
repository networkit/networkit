#ifndef CLIQUE_H_
#define CLIQUE_H_

#include "../graph/Graph.h"
#include "../base/Algorithm.h"
#include <functional>

namespace NetworKit {

class MaximalCliques : public Algorithm {

public:
	MaximalCliques(const Graph& G, bool maximumOnly = false);
	MaximalCliques(const Graph& G, std::function<void(const std::vector<node>&)> callback);

	void run() override;

	const std::vector<std::vector<node>>& getCliques() const;

protected:
	const Graph& G;

	std::vector<std::vector<node>> result;

	std::function<void(const std::vector<node>&)> callback;
	bool maximumOnly;
};

}

#endif /* CLIQUE_H_ */

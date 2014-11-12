/*
 */

#ifndef ADAMICADARATTRIBUTIZER_H_
#define ADAMICADARATTRIBUTIZER_H_

#include "../graph/Graph.h"
#include "AttributeGenerator.h"

namespace NetworKit {

/**
 * An implementation of the triangle counting algorithm by Chiba/Nishizeki.
 */
class AdamicAdarAttributizer : public AttributeGenerator<int, double> {

public:

	virtual std::vector<double> getAttribute(const Graph& graph, const std::vector<int>& attribute) override;

private:
	void removeNode(Graph& graph, node u);
};

} /* namespace NetworKit */

#endif /* ADAMICADARATTRIBUTIZER_H_ */

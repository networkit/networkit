#ifndef TESTSNEMES_H_
#define TESTSNEMES_H_

#include <map>
#include "graph/Graph.h"
#include "centrality/Centrality.h"

namespace NetworKit {

class TestsNemes {

public:
	static double effectiveDiameterExact();
	static double effectiveDiameter(std::string graph, count k);
	static double hopPlot();
	static double closeness();
	static double betweenness();
	static double kpath();
	static double correlation();

private:
	static double calculateCorrelation(const Graph& g, Centrality& centrality1, Centrality& centrality2);

};

} /* namespace NetworKit */

#endif /* TESTSNEMES_H_ */

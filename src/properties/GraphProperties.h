/*
 * GraphProperties.h
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#ifndef GRAPHPROPERTIES_H_
#define GRAPHPROPERTIES_H_


#include "Graph.h"

namespace NetworKit {

class GraphProperties {
public:
	GraphProperties();
	virtual ~GraphProperties();

	static std::vector<count> degreeDistribution(Graph& G);

	static std::vector<double> localClusteringCoefficientPerDegree(Graph& G);

	static std::pair<count, count> minMaxDegree(Graph& G);
};

} /* namespace NetworKit */
#endif /* GRAPHPROPERTIES_H_ */

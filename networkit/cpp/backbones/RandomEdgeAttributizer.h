/*
 * RandomEdgeAttributizer.h
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#ifndef RANDOMEDGEATTRIBUTIZER_H
#define RANDOMEDGEATTRIBUTIZER_H

#include "../edgeproperties/EdgeAttribute.h"

namespace NetworKit {

class RandomEdgeAttributizer : public EdgeAttribute<double> {

public:
	RandomEdgeAttributizer(const Graph& graph, double rneRatio = 0.8);
	
	virtual std::vector<double> getAttribute() override;
	
private:
	const Graph& graph;
	double rneRatio;

};

} /* namespace NetworKit */

#endif // RANDOMEDGEATTRIBUTIZER_H

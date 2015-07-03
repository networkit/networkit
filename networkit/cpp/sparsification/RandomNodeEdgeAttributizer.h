/*
 * RandomNodeEdgeAttributizer.h
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#ifndef RANDOMNODEEDGEATTRIBUTIZER_H
#define RANDOMNODEEDGEATTRIBUTIZER_H

#include "../edgeattributes/EdgeAttribute.h"

namespace NetworKit {

class RandomNodeEdgeAttributizer : public EdgeAttribute<double> {

public:
	RandomNodeEdgeAttributizer(const Graph& graph, double rneRatio = 0.8);

	virtual std::vector<double> getAttribute() override;

private:
	const Graph& graph;
	double rneRatio;

};

} /* namespace NetworKit */

#endif // RANDOMNODEEDGEATTRIBUTIZER_H

/*
 * ChangeCorrectedTriangleAttributizer.h
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#ifndef CHANCECORRECTEDTRIANGLEATTRIBUTIZER_H
#define CHANCECORRECTEDTRIANGLEATTRIBUTIZER_H

#include "../edgeattributes/EdgeAttribute.h"

namespace NetworKit {

class ChanceCorrectedTriangleAttributizer : public EdgeAttribute<double> {

public:
	ChanceCorrectedTriangleAttributizer(const Graph& graph, const std::vector<int>& triangles);
	virtual std::vector<double> getAttribute() override;

private:
	const Graph& graph;
	const std::vector<int>& triangles;

};

}
#endif // CHANCECORRECTEDTRIANGLEATTRIBUTIZER_H

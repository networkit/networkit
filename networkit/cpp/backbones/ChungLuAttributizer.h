/*
 * ChungLuAttributizer.h
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#ifndef CHUNGLUATTRIBUTIZER_H
#define CHUNGLUATTRIBUTIZER_H

#include "../edgeproperties/EdgeAttribute.h"

namespace NetworKit {

class ChungLuAttributizer : public EdgeAttribute<double> {

public:
	ChungLuAttributizer(const Graph& graph);
	virtual std::vector<double> getAttribute() override;

private:
	const Graph& graph;

};

} /* namespace NetworKit */

#endif // CHUNGLUATTRIBUTIZER_H

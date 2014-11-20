/*
 * GeometricMeanAttributizer.h
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#ifndef GEOMETRICMEANATTRIBUTIZER_H
#define GEOMETRICMEANATTRIBUTIZER_H

#include "EdgeAttribute.h"

namespace NetworKit {

class GeometricMeanAttributizer : public EdgeAttribute<double> {

private:
	const Graph& graph;
	const std::vector<double>& attribute;

public:
	GeometricMeanAttributizer(const Graph& graph, const std::vector<double>& attribute);
	virtual std::vector<double> getAttribute() override;
};

} // namespace NetworKit

#endif // GEOMETRICMEANATTRIBUTIZER_H

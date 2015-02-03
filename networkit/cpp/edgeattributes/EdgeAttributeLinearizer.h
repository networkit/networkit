/*
 * EdgeAttributeLinearizer.h
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#ifndef EDGEATTRIBUTELINEARIZER_H
#define EDGEATTRIBUTELINEARIZER_H

#include "EdgeAttribute.h"

namespace NetworKit {

class EdgeAttributeLinearizer : public EdgeAttribute<double> {

private:
	const Graph& graph;
	const std::vector<double>& attribute;
	bool inverse;

public:

	EdgeAttributeLinearizer(const Graph& graph, const std::vector<double>& attribute, bool inverse = false);
	virtual std::vector<double> getAttribute() override;

};

} // namespace NetworKit

#endif // EDGEATTRIBUTELINEARIZER_H

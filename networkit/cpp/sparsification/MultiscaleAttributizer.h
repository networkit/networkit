/*
 * MultiscaleAttributizer.h
 *
 *  Created on: 20.06.2014
 *      Author: Gerd Lindner
 */

#ifndef MULTISCALEATTRIBUTIZER_H_
#define MULTISCALEATTRIBUTIZER_H_

#include "../edgescores/EdgeAttribute.h"

namespace NetworKit {

/**
 * Calculates the multiscale backbone attribute for a given graph. Each edge is
 * assigned the maximum filter value in [0,1] for which the edge will be contained
 * in the multiscale backbone.
 *
 * See "Extracting the multiscale backbone of complex weighted networks" by Serrano et al.
 */
class MultiscaleAttributizer : public EdgeAttribute<double> {

public:

	MultiscaleAttributizer(const Graph& graph, const std::vector<double>& attribute);
	virtual std::vector<double> getAttribute() override;
	double getProbability(count degree, edgeweight normalizedWeight);

private:
	const Graph& graph;
	const std::vector<double>& attribute;

};

}
/* namespace NetworKit */

#endif /* MULTISCALEATTRIBUTIZER_H_ */

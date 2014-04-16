/*
 * AreaWeightedCentrality.h
 *
 *  Created on: 14.04.2014
 *      Author: Maximilian Vogel
 */

#ifndef AREAWEIGHTEDCENTRALITY_H_
#define AREAWEIGHTEDCENTRALITY_H_

#include "Centrality.h"

namespace NetworKit {

/** 
 * Node centrality index which ranks nodes by their degree.
 * Optional normalization by maximum degree.
 */
class AreaWeightedCentrality: public NetworKit::Centrality {
protected:
	std::vector<double> cosine;
public:
	AreaWeightedCentrality(const Graph& G, bool normalized=false);

	void run() override;

	void initCoordinates(std::string lat_file, std::string lon_file);
};

} /* namespace NetworKit */

#endif /* AREAWEIGHTEDCENTRALITY_H_ */

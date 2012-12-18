/*
 * RegionGrowingOverlapper.h
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#ifndef REGIONGROWINGOVERLAPPER_H_
#define REGIONGROWINGOVERLAPPER_H_

#include <vector>

#include "Overlapper.h"
#include "../clustering/Clustering.h"

namespace EnsembleClustering {

class RegionGrowingOverlapper: public EnsembleClustering::Overlapper {

public:

	RegionGrowingOverlapper();

	virtual ~RegionGrowingOverlapper();

	virtual void run(Graph& G, std::vector<Clustering> clusterings);

};

} /* namespace EnsembleClustering */
#endif /* REGIONGROWINGOVERLAPPER_H_ */

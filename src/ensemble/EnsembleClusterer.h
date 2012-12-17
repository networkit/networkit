/*
 * EnsembleClusterer.h
 *
 *  Created on: 17.12.2012
 *      Author: cls
 */

#ifndef ENSEMBLECLUSTERER_H_
#define ENSEMBLECLUSTERER_H_

#include "../clustering/Clusterer.h"

#include <vector>

namespace EnsembleClustering {

class EnsembleClusterer: public EnsembleClustering::Clusterer {

protected:

	std::vector<Clusterer> baseClusterers;

public:

	EnsembleClusterer();

	virtual ~EnsembleClusterer();

	virtual Clustering& run(Graph& G);
};

} /* namespace EnsembleClustering */
#endif /* ENSEMBLECLUSTERER_H_ */

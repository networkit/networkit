/*
 * EnsembleClusterer.h
 *
 *  Created on: 17.12.2012
 *      Author: cls
 */

#ifndef ENSEMBLECLUSTERER_H_
#define ENSEMBLECLUSTERER_H_

#include <vector>


#include "../clustering/Clusterer.h"
#include "../clustering/Clustering.h"
#include "../overlap/Overlapper.h"

namespace EnsembleClustering {

class EnsembleClusterer: public EnsembleClustering::Clusterer {

protected:

	std::vector<Clusterer> baseClusterers;
	std::vector<Clustering> baseClusterings;

public:

	EnsembleClusterer();

	virtual ~EnsembleClusterer();

	virtual Clustering& run(Graph& G);
};

} /* namespace EnsembleClustering */
#endif /* ENSEMBLECLUSTERER_H_ */

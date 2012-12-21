/*
 * EnsembleClusterer.h
 *
 *  Created on: 17.12.2012
 *      Author: cls
 */

#ifndef ENSEMBLECLUSTERER_H_
#define ENSEMBLECLUSTERER_H_

#include <vector>
#include <cmath>


#include "../clustering/Clusterer.h"
#include "../clustering/Clustering.h"
#include "../clustering/Modularity.h"
#include "../clustering/ClusteringGenerator.h"
#include "../overlap/RegionGrowingOverlapper.h"
#include "../coarsening/ClusterContracter.h"

namespace EnsembleClustering {

class EnsembleClusterer: public EnsembleClustering::Clusterer {

protected:

	QualityMeasure* qm;	//!< evaluates clutering quality
	double qBest;	//!< quality of currently best clustering
	Clustering* bestClustering; 	//<! currently best clustering

	Clusterer* finalClusterer;	//!< final clustering algorithm
	std::vector<Clusterer*> baseClusterers;
	std::vector<Clustering*> baseClusterings;

	/**
	 * Check if new clustering is better than previously best clustering.
	 */
	virtual bool isBetterClustering(const Clustering& zeta);

public:

	EnsembleClusterer();

	virtual ~EnsembleClusterer();

	virtual Clustering& run(Graph& G);

	virtual void addBaseClusterer(Clusterer& baseClusterer);

	virtual void setFinalClusterer(Clusterer& clusterer);

};

} /* namespace EnsembleClustering */
#endif /* ENSEMBLECLUSTERER_H_ */

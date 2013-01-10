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
#include <stdexcept>


#include "../clustering/Clusterer.h"
#include "../clustering/Clustering.h"
#include "../clustering/Modularity.h"
#include "../clustering/ClusteringGenerator.h"
#include "../overlap/RegionGrowingOverlapper.h"
#include "../coarsening/ClusterContracter.h"
#include "../coarsening/GraphContraction.h"
#include "../io/GraphIO.h"

namespace EnsembleClustering {

class EnsembleClusterer: public EnsembleClustering::Clusterer {

protected:

	QualityMeasure* qm;	//!< evaluates clutering quality

	Clusterer* finalClusterer;	//!< final clustering algorithm
	std::vector<Clusterer*> baseClusterers;


public:

	EnsembleClusterer();

	virtual ~EnsembleClusterer();

	virtual void addBaseClusterer(Clusterer&  base);

	virtual void setFinalClusterer(Clusterer& final);

	virtual void setQualityMeasure(QualityMeasure& qm);

	virtual Clustering run(Graph& G);


};

} /* namespace EnsembleClustering */
#endif /* ENSEMBLECLUSTERER_H_ */

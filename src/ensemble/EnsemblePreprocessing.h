/*
 * EnsemblePreprocessing.h
 *
 *  Created on: 26.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef ENSEMBLEPREPROCESSING_H_
#define ENSEMBLEPREPROCESSING_H_

#include <vector>
#include <stdexcept>

#include "../clustering/algo/Clusterer.h"

#include "../aux/Log.h"
#include "../Globals.h"
#include "../clustering/base/Clustering.h"
#include "../clustering/base/Modularity.h"
#include "../overlap/HashingOverlapper.h"
#include "../coarsening/ClusterContracter.h"
#include "../coarsening/ClusteringProjector.h"
#include "../clustering/base/RandMeasure.h"
#include "../clustering/base/JaccardMeasure.h"

namespace EnsembleClustering {

class EnsemblePreprocessing: public EnsembleClustering::Clusterer {

protected:

	Clusterer* finalClusterer;	//!< final clustering algorithm
	std::vector<Clusterer*> baseClusterers; //!< ensemble of base clusterers

	Overlapper* overlap; //!< clustering overlap algorithm

public:

	EnsemblePreprocessing();

	virtual ~EnsemblePreprocessing();

	/**
	 * Add a base clusterer to the ensemble.
	 */
	virtual void addBaseClusterer(Clusterer&  base);


	/**
	 * Set final clusterer.
	 */
	virtual void setFinalClusterer(Clusterer& final);

	/**
	 * Set overlap algorithm which combines the results of the base clusterers.
	 */
	virtual void setOverlapper(Overlapper& overlap);

	/**
	 * Run the ensemble clusterer.
	 */
	virtual Clustering run(Graph& G);

	/**
	 * @return string representation.
	 */
	virtual std::string toString() const;

};

} /* namespace EnsembleClustering */
#endif /* ENSEMBLEPREPROCESSING_H_ */

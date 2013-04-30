/*
 * EnsembleMultilevel.h
 *
 *  Created on: 17.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef ENSEMBLEMULTILEVEL_H_
#define ENSEMBLEMULTILEVEL_H_

#include <vector>
#include <cmath>
#include <stdexcept>


#include "../auxiliary/Log.h"
#include "../Globals.h"
#include "../community/Clusterer.h"
#include "../clustering/Clustering.h"
#include "../clustering/Modularity.h"
#include "../clustering/ClusteringGenerator.h"
#include "../overlap/RegionGrowingOverlapper.h"
#include "../overlap/HashingOverlapper.h"
#include "../coarsening/ClusterContracter.h"
#include "../coarsening/GraphContraction.h"
#include "../io/GraphIO.h"
#include "../coarsening/ClusteringProjector.h"
#include "../clustering/NodeStructuralRandMeasure.h"
#include "../clustering/GraphStructuralRandMeasure.h"
#include "../clustering/JaccardMeasure.h"



namespace NetworKit {

class EnsembleMultilevel: public NetworKit::Clusterer {

protected:

	QualityMeasure* qm;	//!< evaluates clutering quality

	Clusterer* finalClusterer;	//!< final clustering algorithm
	std::vector<Clusterer*> baseClusterers; //!< ensemble of base clusterers

	Overlapper* overlap; //!< clustering overlap algorithm

	int h; //!< hierarchy level / iteration counter


public:

	EnsembleMultilevel();

	virtual ~EnsembleMultilevel();

	/**
	 * Add a base clusterer to the ensemble.
	 */
	virtual void addBaseClusterer(Clusterer&  base);

	/**
	 * Set final clusterer.
	 */
	virtual void setFinalClusterer(Clusterer& final);

	/**
	 * Set quality measure to evaluate clusterings.
	 */
	virtual void setQualityMeasure(QualityMeasure& qm);


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

} /* namespace NetworKit */
#endif /* ENSEMBLEMULTILEVEL_H_ */

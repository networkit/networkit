/*
 * DynamicEnsemble.h
 *
 *  Created on: 19.06.2013
 *      Author: cls
 */

#ifndef DYNAMICENSEMBLE_H_
#define DYNAMICENSEMBLE_H_

#include "DynamicCommunityDetector.h"

#include "../community/EnsemblePreprocessing.h"

namespace NetworKit {

class DynamicEnsemble: public NetworKit::DynamicCommunityDetector {
public:

	DynamicEnsemble();

	virtual ~DynamicEnsemble();

	/**
	 * Set the Graph instance. Needs to be called before calling run().
	 */
	virtual void setGraph(Graph& G);

	virtual Clustering run();

	virtual std::string toString() const;

	/** EVENT HANDLING **/

	virtual void onNodeAddition(node u);

	virtual void onNodeRemoval(node u);

	virtual void onEdgeAddition(node u, node v, edgeweight w = 1.0);

	virtual void onEdgeRemoval(node u, node v, edgeweight w = 1.0);

	virtual void onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew);

	virtual void onTimeStep();

	/** ENSEMBLE CONFIGURATION **/

	/**
	 * Add a base clusterer to the ensemble.
	 */
	virtual void addBaseAlgorithm(DynamicCommunityDetector& base);


	/**
	 * Set final clusterer.
	 */
	virtual void setFinalAlgorithm(Clusterer& final);

	/**
	 * Set overlap algorithm which combines the results of the base clusterers.
	 */
	virtual void setOverlapper(Overlapper& overlap);

protected:


	Clusterer* finalAlgo;	//!< final clustering algorithm
	std::vector<DynamicCommunityDetector*> baseAlgos; //!< ensemble of dynamic base algorithms

	Overlapper* overlapAlgo; //!< clustering overlap algorithm

};

} /* namespace NetworKit */
#endif /* DYNAMICENSEMBLE_H_ */

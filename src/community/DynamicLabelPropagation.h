/*
 * DynamicLabelPropagation.h
 *
 *  Created on: 27.03.2013
 *      Author: cls
 */

#ifndef DYNAMICLABELPROPAGATION_H_
#define DYNAMICLABELPROPAGATION_H_

#include <algorithm>

#include "DynamicClusterer.h"
#include "../auxiliary/Timer.h"
#include "PrepStrategy.h"

namespace NetworKit {

typedef cluster label;

class DynamicLabelPropagation: public NetworKit::DynamicClusterer {


public:

	/**
	 * @param[in]	G		graph
	 * @param[in]	theta	update threshold
	 */
	DynamicLabelPropagation(Graph& G, count theta, std::string strategy = "reactivate");

	virtual ~DynamicLabelPropagation();

	virtual Clustering run();

	virtual std::string toString() const;

	// GRAPH DYNAMICS INTERFACE

	virtual void onNodeAddition(node u);

	virtual void onNodeRemoval(node u);

	virtual void onEdgeAddition(node u, node v);

	virtual void onEdgeRemoval(node u, node v);

	virtual void onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew);

protected:

	Clustering labels;					//!< the labelling/clustering
	std::vector<bool> activeNodes;		//!< which nodes are currently active?
	std::vector<double> weightedDegree; //!< precompute and update weighted degree for performance reasons
	count updateThreshold;
	count nUpdated; 					//!< number of nodes updated in last iteration (?)

	PrepStrategy* prepStrategy;			//!< a prep strategy reacts to graph events by preparing the algorithm's clustering


	// PREP STRATEGIES
	// prep strategies encapsulate different update schemes and make DynamicLabelPropagation modular

	class Reactivate : public NetworKit::PrepStrategy {

		friend class DynamicLabelPropagation;

	public:

		Reactivate(DynamicLabelPropagation* dynPLP);

		virtual ~Reactivate();

		/**
		 * New node u becomes a singleton and active.
		 */
		virtual void onNodeAddition(node u);

		/**
		 * Removed node u becomes permanently inactive.
		 */
		virtual void onNodeRemoval(node u);

		/**
		 * Reactivate u and v on addition of edge {u,v}
		 */
		virtual void onEdgeAddition(node u, node v);

		/**
		 * Reactivate u and v on removal of edge {u,v}
		 */
		virtual void onEdgeRemoval(node u, node v);

		/**
		 * Same reaction as onEdgeAddition
		 */
		virtual void onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew);

	protected:

		DynamicLabelPropagation* dynPLP;

	};


	class ReactivateNeighbors : public NetworKit::PrepStrategy {

		friend class DynamicLabelPropagation;

	public:

		ReactivateNeighbors(DynamicLabelPropagation* dynPLP);

		virtual ~ReactivateNeighbors();

		/**
		 * New node u becomes a singleton and active.
		 */
		virtual void onNodeAddition(node u);

		/**
		 * Removed node u becomes permanently inactive.
		 */
		virtual void onNodeRemoval(node u);

		/**
		 * Reactivate u and v on addition of edge {u,v}
		 */
		virtual void onEdgeAddition(node u, node v);

		/**
		 * Reactivate u and v on removal of edge {u,v}
		 */
		virtual void onEdgeRemoval(node u, node v);

		/**
		 * Same reaction as onEdgeAddition
		 */
		virtual void onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew);

	protected:

		DynamicLabelPropagation* dynPLP;

	};

};

} /* namespace NetworKit */
#endif /* DYNAMICLABELPROPAGATION_H_ */

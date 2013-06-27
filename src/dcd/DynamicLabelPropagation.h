/*
 * DynamicLabelPropagation.h
 *
 *  Created on: 27.03.2013
 *      Author: cls
 */

#ifndef DYNAMICLABELPROPAGATION_H_
#define DYNAMICLABELPROPAGATION_H_

#include <algorithm>

#include "DynamicCommunityDetector.h"
#include "../auxiliary/Timer.h"
#include "PrepStrategy.h"

namespace NetworKit {

typedef cluster label;

class DynamicLabelPropagation: public NetworKit::DynamicCommunityDetector {


public:

	DynamicLabelPropagation(); // nullary constructor needed for Python interface - do not use this to construct instance

	/**
	 * @param[in]	theta	update threshold
	 */
	DynamicLabelPropagation(count theta = 0, std::string strategy = "Reactivate");

	virtual ~DynamicLabelPropagation();


	/**
	 * Set the input graph and initialize data structures.
	 */
	virtual void setGraph(Graph& G);

	/**
	 * Run the Label Propagation community detection algorithm and produce a clustering.
	 * This reuses the previous clustering and takes collected graph modifications
	 * into account.
	 */
	virtual Clustering run();

	virtual std::string toString() const;

	// GRAPH DYNAMICS INTERFACE

	virtual void onNodeAddition(node u);

	virtual void onNodeRemoval(node u);

	virtual void onEdgeAddition(node u, node v, edgeweight w = 1.0);

	virtual void onEdgeRemoval(node u, node v, edgeweight w = 1.0);

	virtual void onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew);

	virtual void onTimeStep();


protected:

	count updateThreshold;
	Clustering labels;					//!< the labelling/clustering
	std::vector<bool> activeNodes;		//!< which nodes are currently active?
	std::vector<double> weightedDegree; //!< precompute and update weighted degree for performance reasons
	count nUpdated; 					//!< number of nodes updated in last iteration (?)

	PrepStrategy* prepStrategy;			//!< a prep strategy reacts to graph events by preparing the algorithm's clustering

	// main timer is in superclass
//	Aux::Timer prepTimer; 				//!< measure time spent in prep strategy
//	std::vector<count> prepTimerHistory;//!< record time spent in prep strategy

	// PREP STRATEGIES
	// prep strategies encapsulate different update schemes and make DynamicLabelPropagation modular

	/**
	 * Reactivates nodes affected by a change.
	 */
	class Reactivate : public NetworKit::PrepStrategy {

		friend class DynamicLabelPropagation;

	public:

		Reactivate(DynamicLabelPropagation* dynPLP);

		virtual ~Reactivate();

		std::string toString();

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
		virtual void onEdgeAddition(node u, node v, edgeweight w = 1.0);

		/**
		 * Reactivate u and v on removal of edge {u,v}
		 */
		virtual void onEdgeRemoval(node u, node v, edgeweight w = 1.0);

		/**
		 * Same reaction as onEdgeAddition
		 */
		virtual void onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew);

		/**
		 * Ignore time steps.
		 */
		virtual void onTimeStep();

	protected:

		DynamicLabelPropagation* dynPLP;

	};


	/**
	 * Activates nodes affected by a change and also their neighbors.
	 */
	class ReactivateNeighbors : public NetworKit::PrepStrategy {

		friend class DynamicLabelPropagation;

	public:

		ReactivateNeighbors(DynamicLabelPropagation* dynPLP);

		virtual ~ReactivateNeighbors();

		virtual std::string toString();

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
		virtual void onEdgeAddition(node u, node v, edgeweight w = 1.0);

		/**
		 * Reactivate u and v on removal of edge {u,v}
		 */
		virtual void onEdgeRemoval(node u, node v, edgeweight w = 1.0);

		/**
		 * Same reaction as onEdgeAddition
		 */
		virtual void onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew);

		/**
		 * Ignore time step.
		 */
		virtual void onTimeStep();

	protected:

		DynamicLabelPropagation* dynPLP;

	};


	/**
	 * Turn nodes affected by a change into singletons. Since this is a label change,
	 * nodes in the neighborhood get activated.
	 */
	class Isolate : public NetworKit::PrepStrategy {

			friend class DynamicLabelPropagation;

		public:

			Isolate(DynamicLabelPropagation* dynPLP);

			virtual ~Isolate();

			std::string toString();

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
			virtual void onEdgeAddition(node u, node v, edgeweight w = 1.0);

			/**
			 * Reactivate u and v on removal of edge {u,v}
			 */
			virtual void onEdgeRemoval(node u, node v, edgeweight w = 1.0);

			/**
			 * Same reaction as onEdgeAddition
			 */
			virtual void onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew);

			/**
			 * Ignore time steps.
			 */
			virtual void onTimeStep();

		protected:

			DynamicLabelPropagation* dynPLP;

		};


	/**
	 * Turn nodes affected by a change and their neighbors into singletons. Since this is a label change,
	 * nodes in the 2-neighborhood get activated.
	 */
	class IsolateNeighbors : public NetworKit::PrepStrategy {

			friend class DynamicLabelPropagation;

		public:

			IsolateNeighbors(DynamicLabelPropagation* dynPLP);

			virtual ~IsolateNeighbors();

			std::string toString();

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
			virtual void onEdgeAddition(node u, node v, edgeweight w = 1.0);

			/**
			 * Reactivate u and v on removal of edge {u,v}
			 */
			virtual void onEdgeRemoval(node u, node v, edgeweight w = 1.0);

			/**
			 * Same reaction as onEdgeAddition
			 */
			virtual void onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew);

			/**
			 * Ignore time steps.
			 */
			virtual void onTimeStep();

		protected:

			DynamicLabelPropagation* dynPLP;

		};

};

} /* namespace NetworKit */
#endif /* DYNAMICLABELPROPAGATION_H_ */

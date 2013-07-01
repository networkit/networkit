/*
 * TDynamicLabelPropagation.h
 *
 *  Created on: 21.06.2013
 *      Author: cls
 */

#ifndef TDYNAMICLABELPROPAGATION_H_
#define TDYNAMICLABELPROPAGATION_H_

#include "DynamicCommunityDetector.h"

namespace NetworKit {

typedef cluster label;

template <class PrepStrategy> class TDynamicLabelPropagation: public DynamicCommunityDetector {

	friend class Isolate; // PrepStrategies must be friend classes to access internal data
	friend class IsolateNeighbors;

public:
	TDynamicLabelPropagation(count theta=0);

	virtual ~TDynamicLabelPropagation();

	/**
	 * Set the Graph instance. Needs to be called before calling run().
	 */
	virtual void setGraph(Graph& G) override;

	virtual Clustering run() override;

	virtual std::string toString() const override;

	virtual void onNodeAddition(node u) override;

	virtual void onNodeRemoval(node u) override;

	virtual void onEdgeAddition(node u, node v, edgeweight w = 1.0) override;

	virtual void onEdgeRemoval(node u, node v, edgeweight w = 1.0) override;

	virtual void onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew) override;

	virtual void onTimeStep() override;


protected:

	PrepStrategy prepStrategy; 					//!< instance of the prep strategy

	count updateThreshold;
	Clustering labels;					//!< the labelling/clustering
	std::vector<bool> activeNodes;		//!< which nodes are currently active?
	std::vector<double> weightedDegree; //!< precompute and update weighted degree for performance reasons
	count nUpdated; 					//!< number of nodes updated in last iteration (?)
};




/**
 * Turn nodes affected by a change into singletons. Since this is a label change,
 * nodes in the neighborhood get activated.
 */
class Isolate {

public:
	Isolate() : dynPLP(NULL) {

	};

	// TODO: not necessary, can be done in the constructor
	void setAlgorithm(TDynamicLabelPropagation<Isolate>* instance) {
		this->dynPLP = instance;
	}

	inline void onNodeAddition(node u) {
		// dynPLP has already made the new node a singleton
		dynPLP->G->forNeighborsOf(u, [&](node v) {
			dynPLP->activeNodes[v] = true;
		});
	};

	inline void onNodeRemoval(node u) {
		assert (dynPLP->activeNodes[u] == false); // dynPLP has permanently deactivated the node
		// incident edges must have been removed before, so neighborhood is already active at this point
	};

	inline void onEdgeAddition(node u, node v, edgeweight w = 1.0) {
		dynPLP->labels.toSingleton(u);
		dynPLP->activeNodes[u] = true;
		dynPLP->labels.toSingleton(v);
		dynPLP->activeNodes[v] = true;
		dynPLP->G->forNeighborsOf(u, [&](node x) {
			dynPLP->activeNodes[x] = true;
		});
		dynPLP->G->forNeighborsOf(v, [&](node x) {
			dynPLP->activeNodes[x] = true;
		});

	};

	inline void onEdgeRemoval(node u, node v, edgeweight w = 1.0) {
		this->onEdgeAddition(u, v, w);
	};

	inline void onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew) {
		this->onEdgeAddition(u, v, wNew);
	};

	inline void onTimeStep() {
		// do nothing
	};

	inline std::string toString() const {
		return "Isolate";
	}

protected:

	TDynamicLabelPropagation<Isolate>* dynPLP;
};



/**
 * Turn nodes affected by a change into singletons. Since this is a label change,
 * nodes in the neighborhood get activated.
 */
class IsolateNeighbors {

public:
	IsolateNeighbors() : dynPLP(NULL) {

	};

	// TODO: not necessary, can be done in the constructor
	void setAlgorithm(TDynamicLabelPropagation<IsolateNeighbors>* instance) {
		this->dynPLP = instance;
	}

	inline void onNodeAddition(node u) {
		assert (dynPLP->activeNodes[u] == true);
		assert (dynPLP->G->degree(u) == 0); // new node has no incident edges
	};

	inline void onNodeRemoval(node u) {
		assert (dynPLP->activeNodes[u] == false); // dynPLP has permanently deactivated the node
		// incident edges must have been removed before, so neighborhood is already active at this point
	};

	inline void onEdgeAddition(node u, node v, edgeweight w = 1.0) {
		// turn u, v and their neighbors into singletons, activate 2-neighborhood
		dynPLP->labels.toSingleton(u);
		dynPLP->activeNodes[u] = true;
		dynPLP->labels.toSingleton(v);
		dynPLP->activeNodes[v] = true;
		dynPLP->G->forNeighborsOf(u, [&](node x) {
			dynPLP->labels.toSingleton(x);
			dynPLP->G->forNeighborsOf(x, [&](node y){
				dynPLP->activeNodes[y] = true;
			});
		});
		dynPLP->G->forNeighborsOf(v, [&](node x) {
			dynPLP->labels.toSingleton(x);
			dynPLP->G->forNeighborsOf(x, [&](node y){
				dynPLP->activeNodes[y] = true;
			});
		});

	};

	inline void onEdgeRemoval(node u, node v, edgeweight w = 1.0) {
		this->onEdgeAddition(u, v, w);
	};

	inline void onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew) {
		this->onEdgeAddition(u, v, wNew);
	};

	inline void onTimeStep() {
		// do nothing
	};

	inline std::string toString() const {
		return "IsolateNeighbors";
	}

protected:

	TDynamicLabelPropagation<IsolateNeighbors>* dynPLP;
};

template<class PrepStrategy>
inline TDynamicLabelPropagation<PrepStrategy>::TDynamicLabelPropagation(count theta) : DynamicCommunityDetector(), updateThreshold(theta), nUpdated(0) {
	this->prepStrategy.setAlgorithm(this);
}

template<class PrepStrategy>
inline TDynamicLabelPropagation<PrepStrategy>::~TDynamicLabelPropagation() {
}

template<class PrepStrategy>
inline void TDynamicLabelPropagation<PrepStrategy>::setGraph(Graph& G) {
	if (!G.isEmpty()) {
		throw std::runtime_error("G is not an empty graph. Currently, it is assumed that this algorithm is initialized with an empty graph, which is then constructed incrementally");
	}

	this->G = &G;
}

template<class PrepStrategy>
inline Clustering TDynamicLabelPropagation<PrepStrategy>::run() {
	if (this->G == NULL) {
		throw std::runtime_error("pointer to current graph was not initialized - call setGraph first");
	}

	INFO("running TDynamicLabelPropagation at t=" << G->time());

	count nIterations = 0;
	nUpdated = G->numberOfNodes(); // starts while loop - TODO: correct?

	runtime.start();
	while (nUpdated > updateThreshold) {
		nUpdated = 0;
		nIterations += 1;
		G->parallelForNodes([&](node u){
			if ((activeNodes[u]) && (G->degree(u) > 0)) {
				std::map<label, double> labelWeights; // neighborLabelCounts maps label -> frequency in the neighbors

				// weigh the labels in the neighborhood of v
				G->forWeightedNeighborsOf(u, [&](node v, edgeweight weight) {
					label lv = labels[v];
					labelWeights[lv] += weight; // add weight of edge {v, w}
				});

				// get heaviest label
				label heaviest = std::max_element(labelWeights.begin(),
								labelWeights.end(),
								[](const std::pair<label, edgeweight>& p1, const std::pair<label, edgeweight>& p2) {
									return p1.second < p2.second;})->first;

				if (labels[u] != heaviest) { // UPDATE
					labels[u] = heaviest;
					nUpdated += 1; // TODO: atomic update?
					G->forNeighborsOf(u, [&](node v) {
						activeNodes[v] = true;
					});
				} else {
					activeNodes[u] = false;
				}
			} else { /* node is isolated */ }
		}); // end parallel for nodes
		TRACE("nUpdated = " << nUpdated);
	} // end while

	runtime.stop();
	this->timerHistory.push_back(runtime.elapsed().count());
	INFO("[DONE] iteration #" << nIterations << ", time spent: " << runtime.elapsedTag());


	return labels;
}

template<class PrepStrategy>
inline std::string TDynamicLabelPropagation<PrepStrategy>::toString() const {
	std::stringstream strm;
	strm << "TDynamicLabelPropagation(" << this->prepStrategy.toString() << ")";
	return strm.str();
}

template<class PrepStrategy>
inline void TDynamicLabelPropagation<PrepStrategy>::onNodeAddition(node u) {
	// update data structures
	activeNodes.push_back(true); // new node is active
	weightedDegree.push_back(0.0);
	labels.append(u); // extend label array by 1 entry
	labels.toSingleton(u); // create singleton

	prepStrategy.onNodeAddition(u);
}

template<class PrepStrategy>
inline void TDynamicLabelPropagation<PrepStrategy>::onNodeRemoval(node u) {
	assert (G->degree(u) == 0);
	assert (weightedDegree[u] == 0.0);

	this->activeNodes[u] = false; // assumption: this node can never be reactivated

	prepStrategy.onNodeRemoval(u);
}

template<class PrepStrategy>
inline void TDynamicLabelPropagation<PrepStrategy>::onEdgeAddition(
		node u, node v, edgeweight w) {
	// update weighted degree
	if (u != v) {
		weightedDegree[u] += w;
		weightedDegree[v] += w;
	} else {
		weightedDegree[u] += w; // self-loop case
	}

	// assert that this is consistent with the graph
	assert (G->weightedDegree(u) == weightedDegree[u]);
	assert (G->weightedDegree(v) == weightedDegree[v]);
	// apply prep strategy
	this->prepStrategy.onEdgeAddition(u, v);
}

template<class PrepStrategy>
inline void TDynamicLabelPropagation<PrepStrategy>::onEdgeRemoval(
		node u, node v, edgeweight w) {
	// update weighted degree
	if (u != v) {
		weightedDegree[u] -= w;
		weightedDegree[v] -= w;
	} else {
		weightedDegree[u] -= w; // self-loop case
	}
	// assert that this is consistent with the graph
	assert (G->weightedDegree(u) == weightedDegree[u]);
	assert (G->weightedDegree(v) == weightedDegree[v]);
	// apply prep strategy
	this->prepStrategy.onEdgeRemoval(u, v);
}

template<class PrepStrategy>
inline void TDynamicLabelPropagation<PrepStrategy>::onWeightUpdate(
		node u, node v, edgeweight wOld, edgeweight wNew) {
	//update weighted degree
	if (u != v) {
		weightedDegree[u] += (wNew - wOld);
		weightedDegree[v] += (wNew - wOld);
	} else {
		weightedDegree[u] += (wNew - wOld);
	}

	this->prepStrategy.onWeightUpdate(u, v, wOld, wNew);
}

template<class PrepStrategy>
inline void TDynamicLabelPropagation<PrepStrategy>::onTimeStep() {
	// do nothing
}


} /* namespace NetworKit */


#endif /* TDYNAMICLABELPROPAGATION_H_ */

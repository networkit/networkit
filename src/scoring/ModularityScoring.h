/*
 * Modularity.h
 *
 *  Created on: 15.10.2012
 *      Author: cls
 */

#ifndef MODULARITY_H_
#define MODULARITY_H_

#include "EdgeScoring.h"

#include "../clustering/base/Clustering.h"


namespace EnsembleClustering {

// TODO: implement modularity as in Python prototype

template<typename T>
class ModularityScoring: public EnsembleClustering::EdgeScoring<T> {

protected:

	double totalEdgeWeight;	//!< total weight of the graph

public:

	/**
	 * @param[in]	G	a graph instance
	 *
	 * Do not modify the graph while using this instance of ModularityScoring.
	 */
	ModularityScoring(Graph& G);

	virtual ~ModularityScoring();


	virtual void scoreEdges(int attrId);


	/**
	 * Returns an edge score for an edge (u,v) which expresses the
	 * modularity increase which can be gained by merging
	 * the clusters of u and v.
	 *
	 *		 $$\Delta mod(c, d) := \frac{1}{2 \omega(E)} \left ( 2 \omega(E) \omega(c,d) - \omega(c) \omega(d) \right ) $$
	 *
	 * @param[in]	u	source node id
	 * @param[out]	v	target node id
	 *
	 */
	virtual T edgeScore(node u, node v) const;



//	/**
//	 * Calculates the difference in modularity that would result from a merger of
//	 * two clusters.
//	 *
//	 */
//	virtual double deltaMod(cluster c, cluster d) =0;
//
//	virtual double cutweight(cluster c, cluster d) =0;
//
//	virtual double weight(cluster c) =0;
};


template<typename T>
ModularityScoring<T>::ModularityScoring(Graph& G) : EdgeScoring<T>(G) {
	this->totalEdgeWeight = this->G->totalEdgeWeight();
}

template<typename T>
ModularityScoring<T>::~ModularityScoring() {
	// TODO Auto-generated destructor stub
}

template<typename T>
T ModularityScoring<T>::edgeScore(node u, node v) const {
	// return G->
	// TODO
	// calculate $$\Delta mod(c, d) := \frac{1}{2 \omega(E)} \left ( 2 \omega(E) \omega(c,d) - \omega(c) \omega(d) \right ) $$

	double deltaMod = 0.0; // (1 / (2 * omegaE)) * (2 * omegaE * G->weight(u, v) - G->weight(u) * G->weight(v));
	return deltaMod;
}

template<typename T>
void EnsembleClustering::ModularityScoring<T>::scoreEdges(int attrId) {
	this->G->forEdgesWithAttribute_double(attrId, [&](node u, node v, double attr) {
		attr = edgeScore(u, v);
		this->G->setAttribute_double(u, v, attrId, attr);
	});
}


} /* namespace EnsembleClustering */

#endif /* MODULARITY_H_ */

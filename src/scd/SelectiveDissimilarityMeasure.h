/*
 * SelectiveDissimilarityMeasure.h

 *
 *  Created on: 18.06.2013
 *      Author: cls, Yassine Marrakchi
 */

#ifndef SELECTIVEDISSIMILARITYMEASURE_H_
#define SELECTIVEDISSIMILARITYMEASURE_H_

#include <unordered_set>
#include <unordered_map>

#include "../clustering/Clustering.h"
#include "../clustering/Clustering.h"


namespace NetworKit {

/**
 * Abstract base class for all measures which evaluate a single community
 * with respect to a given ground-truth clustering of the graph.
 */
class SelectiveDissimilarityMeasure {

public:

	SelectiveDissimilarityMeasure();

	virtual ~SelectiveDissimilarityMeasure();

	/**
	 * @param[in] seedNode		seed node
	 * @param[in] community		pointer to the considered community
	 * @param[in] groundTruth 	pointer to the ground truth clustering
	 *
	 * return dissimilarity value
	 */
	virtual double localDissimilarity(const node seedNode,const std::unordered_set<node>& community, const Clustering& groundTruth) = 0;

};

/**
 * Jaccard index as selective similarity measure
 */
class JaccardIndex : SelectiveDissimilarityMeasure {

public:

	JaccardIndex();

	virtual ~JaccardIndex();

	virtual double localDissimilarity(const node seedNode,const std::unordered_set<node>& community, const Clustering& groundTruth);

};

/**
 * Precision as selective similarity measure
 */
class Precision : SelectiveDissimilarityMeasure {

public:

	Precision();

	virtual ~Precision();

	double localDissimilarity(const node seedNode,const std::unordered_set<node>& community, const Clustering& groundTruth);

};

/**
 * Recall as selective similarity measure
 */
class Recall : SelectiveDissimilarityMeasure {

public:

	Recall();

	virtual ~Recall();

	double localDissimilarity(const node seedNode,const std::unordered_set<node>& community, const Clustering& groundTruth);

};

/**
 * NMI as selective similarity measure
 */
class NMI : SelectiveDissimilarityMeasure {

public:

	NMI();

	virtual ~NMI();

	double localDissimilarity(const node seedNode,const std::unordered_set<node>& community, const Clustering& groundTruth);

};




} /* namespace NetworKit */
#endif /* SELECTIVEDISSIMILARITYMEASURE_H_ */

/*
 * SelectiveDissimilarityMeasure.h
 *
 *  Created on: 18.06.2013
 *      Author: cls
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

	virtual double localDissimilarity(const node seedNode,const std::unordered_set<node>& community, const Clustering& groundTruth) = 0;

	virtual double getDissimilarity(const std::unordered_map<node, std::unordered_set<node> > communities, const Clustering& groundTruth) = 0;
};

class JaccardIndex : SelectiveDissimilarityMeasure {

public:

	JaccardIndex();

	virtual ~JaccardIndex();

	virtual double localDissimilarity(const node seedNode,const std::unordered_set<node>& community, const Clustering& groundTruth);

	virtual double getDissimilarity(const std::unordered_map<node, std::unordered_set<node> > communities, const Clustering& groundTruth);
};

class Precision : SelectiveDissimilarityMeasure {

public:

	Precision();

	virtual ~Precision();

	double localDissimilarity(const node seedNode,const std::unordered_set<node>& community, const Clustering& groundTruth);

	virtual double getDissimilarity(const std::unordered_map<node, std::unordered_set<node> > communities, const Clustering& groundTruth);
};

class Recall : SelectiveDissimilarityMeasure {

public:

	Recall();

	virtual ~Recall();

	double localDissimilarity(const node seedNode,const std::unordered_set<node>& community, const Clustering& groundTruth);

	virtual double getDissimilarity(const std::unordered_map<node, std::unordered_set<node> > communities, const Clustering& groundTruth);
};

class NMI : SelectiveDissimilarityMeasure {

public:

	NMI();

	virtual ~NMI();

	double localDissimilarity(const node seedNode,const std::unordered_set<node>& community, const Clustering& groundTruth);

	virtual double getDissimilarity(const std::unordered_map<node, std::unordered_set<node> > communities, const Clustering& groundTruth);

};




} /* namespace NetworKit */
#endif /* SELECTIVEDISSIMILARITYMEASURE_H_ */

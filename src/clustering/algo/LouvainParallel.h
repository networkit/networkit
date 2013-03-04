/*
 * LouvainParallel.h
 *
 *  Created on: 27.02.2013
 *      Author: cls
 */

#ifndef LOUVAINPARALLEL_H_
#define LOUVAINPARALLEL_H_

#include "Clusterer.h"

#include "../../coarsening/ClusterContracter.h"
#include "../../coarsening/ClusteringProjector.h"
#include "../../base/IndexMap.h"

namespace EnsembleClustering {

class LouvainParallel: public EnsembleClustering::Clusterer {

protected:
	bool anyChange;	//!< indicates whether any change was made to the clustering in the last pass over the nodes

public:

	LouvainParallel();

	virtual ~LouvainParallel();

	virtual Clustering pass(Graph& G);

	virtual Clustering run(Graph& G);

	/**
	 * @return string representation of algorithm and parameters.
	 */
	virtual std::string toString() const;
};

} /* namespace EnsembleClustering */
#endif /* LOUVAINPARALLEL_H_ */

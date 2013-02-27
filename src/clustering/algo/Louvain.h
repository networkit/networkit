/*
 * Louvain.h
 *
 *  Created on: 25.02.2013
 *      Author: cls
 */

#ifndef LOUVAIN_H_
#define LOUVAIN_H_

#include "Clusterer.h"
#include "../../coarsening/ClusterContracter.h"
#include "../../coarsening/ClusteringProjector.h"
#include "../../base/IndexMap.h"

namespace EnsembleClustering {

class Louvain: public EnsembleClustering::Clusterer {

protected:
	bool anyChange;	//!< indicates whether any change was made to the clustering in the last pass over the nodes


public:
	Louvain();

	virtual ~Louvain();

	virtual Clustering pass(Graph& G);

	virtual Clustering run(Graph& G);

	/**
	 * @return string representation of algorithm and parameters.
	 */
	virtual std::string toString() const;
};

} /* namespace EnsembleClustering */
#endif /* LOUVAIN_H_ */

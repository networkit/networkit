/*
 * Louvain.h
 *
 *  Created on: 25.02.2013
 *      Author: cls
 */

#ifndef LOUVAIN_H_
#define LOUVAIN_H_

#include "Clusterer.h"

namespace EnsembleClustering {

class Louvain: public EnsembleClustering::Clusterer {
public:
	Louvain();
	virtual ~Louvain();

	virtual Clustering run(Graph& G);
};

} /* namespace EnsembleClustering */
#endif /* LOUVAIN_H_ */

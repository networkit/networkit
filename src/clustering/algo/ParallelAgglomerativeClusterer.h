/*
 * ParallelAgglomerativeClusterer.h
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#ifndef PARALLELAGGLOMERATIVECLUSTERER_H_
#define PARALLELAGGLOMERATIVECLUSTERER_H_

#include "Clusterer.h"

namespace EnsembleClustering {

class ParallelAgglomerativeClusterer: public EnsembleClustering::Clusterer {

public:

	ParallelAgglomerativeClusterer();

	virtual ~ParallelAgglomerativeClusterer();
};

} /* namespace EnsembleClustering */
#endif /* SCOREMATCHCONTRACT_H_ */

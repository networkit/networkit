/*
 * RandMeasure.h
 *
 *  Created on: 19.01.2013
 *      Author: cls
 */

#ifndef RANDMEASURE_H_
#define RANDMEASURE_H_

#include "DissimilarityMeasure.h"

namespace EnsembleClustering {

class RandMeasure: public EnsembleClustering::DissimilarityMeasure {

public:

	RandMeasure();

	virtual ~RandMeasure();

	virtual double getDissimilarity(Graph& G, Clustering& first, Clustering& second);

};

} /* namespace EnsembleClustering */
#endif /* RANDMEASURE_H_ */

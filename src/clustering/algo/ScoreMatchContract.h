/*
 * ScoreMatchContract.h
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#ifndef SCOREMATCHCONTRACT_H_
#define SCOREMATCHCONTRACT_H_

#include "Clusterer.h"

namespace EnsembleClustering {

class ScoreMatchContract: public EnsembleClustering::Clusterer {
public:
	ScoreMatchContract();
	virtual ~ScoreMatchContract();
};

} /* namespace EnsembleClustering */
#endif /* SCOREMATCHCONTRACT_H_ */

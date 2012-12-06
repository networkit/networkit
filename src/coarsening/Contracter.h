/*
 * Contracter.h
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#ifndef CONTRACTER_H_
#define CONTRACTER_H_


namespace EnsembleClustering {

#include "../graph/Graph.h"

class Contracter {

public:

	Contracter();

	virtual ~Contracter();

	virtual node contract(node u, node v);
};


} // namespace


#endif /* CONTRACTER_H_ */

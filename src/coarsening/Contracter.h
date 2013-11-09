/*
 * Contracter.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef CONTRACTER_H_
#define CONTRACTER_H_

namespace NetworKit {

/**
 * Abstract base class for graph coarsening/contraction algorithms.
 */
class Contracter {

public:

	Contracter();

	virtual ~Contracter();


};


} // namespace


#endif /* CONTRACTER_H_ */

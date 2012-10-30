/*
 * Contracter.h
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#ifndef CONTRACTER_H_
#define CONTRACTER_H_


namespace EnsembleClustering {

class Contracter {

public:

	Contracter();

	virtual ~Contracter();

	virtual Node contract(Node u, Node v);
};


} // namespace


#endif /* CONTRACTER_H_ */

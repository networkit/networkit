/*
* DorogovtsevMendesGenerator.h
*
*  Created on: 27.05.2014
*      Author: Christian Staudt
*/

#ifndef DOROGOVTSEVMENDESGENERATOR_H_
#define DOROGOVTSEVMENDESGENERATOR_H_

#include "StaticGraphGenerator.h"

namespace NetworKit {

/**
 * @ingroup generators
 */
class DorogovtsevMendesGenerator: public StaticGraphGenerator {

public:
	/**
	* TODO:
	*
	* @param nNodes 	number of nodes in target graph
	*/
	DorogovtsevMendesGenerator(count nNodes);

	virtual Graph generate();

protected:
		count nNodes;

};

} /* namespace NetworKit */
#endif /* DOROGOVTSEVMENDESGENERATOR_H_ */

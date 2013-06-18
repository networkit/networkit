/*
 * PseudoDynamic.h
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#ifndef PSEUDODYNAMIC_H_
#define PSEUDODYNAMIC_H_

#include "../generators/DynamicGraphSource.h"

namespace NetworKit {

class PseudoDynamic: public NetworKit::DynamicGraphSource {
public:

	PseudoDynamic(const Graph& G);

	virtual ~PseudoDynamic();

	virtual void initializeGraph();

	virtual void generate();


protected:

	const Graph& G;
};

} /* namespace NetworKit */
#endif /* PSEUDODYNAMIC_H_ */

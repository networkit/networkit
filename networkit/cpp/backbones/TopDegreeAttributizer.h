/*
 * TopDegreeAttributizer.h
 *
 *  Created on: 26.08.2014
 *      Author: Gerd Lindner
 */

#ifndef TOPDEGREEATTRIBUTIZER_H_
#define TOPDEGREEATTRIBUTIZER_H_

#include "BackboneCalculator.h"

namespace NetworKit {

// TODO: documentation
/**
 * Experimental
 */
class TopDegreeAttributizer : public AttributeGenerator<int, count> {

public:

	TopDegreeAttributizer(); // TODO: declaring default constructors is unneccessary, right?
	~TopDegreeAttributizer() = default; // TODO: declaring destructors is mostly unnecessary
	std::vector<count> getAttribute(const Graph& graph, const std::vector<int>& attribute);

};

}
/* namespace NetworKit */
#endif /* TOPDEGREEATTRIBUTIZER_H_ */

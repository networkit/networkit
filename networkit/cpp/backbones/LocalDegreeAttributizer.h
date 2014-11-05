/*
 * LocalDegreeAttributizer.h
 *
 *  Created on: 28.08.2014
 *      Author: Gerd Lindner
 */

#ifndef LOCALDEGREEATTRIBUTIZER_H_
#define LOCALDEGREEATTRIBUTIZER_H_

#include "AttributeGenerator.h"

namespace NetworKit {

/**
 * EXPERIMENTAL
 */
class LocalDegreeAttributizer : public AttributeGenerator<int, double> {

public:

	LocalDegreeAttributizer();

	std::vector<double> getAttribute(const Graph& graph, const std::vector<int>& attribute);

};

}
/* namespace NetworKit */
#endif /* LOCALDEGREEATTRIBUTIZER_H_ */

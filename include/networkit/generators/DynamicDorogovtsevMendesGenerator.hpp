/*
 * DynamicDorogovtsevMendesGenerator.h
 *
 *  Created on: 03.02.2014
 *      Author: cls
 */

#ifndef DYNAMICDOROGOVTSEVMENDESGENERATOR_H_
#define DYNAMICDOROGOVTSEVMENDESGENERATOR_H_

#include <networkit/generators/DynamicGraphGenerator.hpp>
#include <networkit/auxiliary/Random.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 */
class DynamicDorogovtsevMendesGenerator: public DynamicGraphGenerator {

public:

	DynamicDorogovtsevMendesGenerator();

	std::vector<GraphEvent> generate(count nSteps) override;

private:

	std::vector<std::pair<node, node> > edges;
	bool initial;
	node u; // current node

};

} /* namespace NetworKit */

#endif /* DYNAMICDOROGOVTSEVMENDESGENERATOR_H_ */

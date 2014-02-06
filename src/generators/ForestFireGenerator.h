/*
 * ForestFireGenerator.h
 *
 *  Created on: 17.01.2014
 *      Author: cls
 */

#ifndef FORESTFIREGENERATOR_H_
#define FORESTFIREGENERATOR_H_

#include "DynamicGraphGenerator.h"

namespace NetworKit {

/**
 * The Forest Fire generative model produces dynamic graphs with the following properties:
 * - heavy tailed degree distribution
 * - communities
 * - densification power law
 * - shringking diameter
 *
 * This class implements a variant of the original model which produces undirected graphs.
 *
 * see Leskovec, Kleinberg, Faloutsos: Graphs over Tim: Densification Laws,
 * 	Shringking Diameters and Possible Explanations
 */
class ForestFireGenerator: public NetworKit::DynamicGraphGenerator {

public:

	/**
	 * @param	p	"burning probability", controls the amount of connections formed by each new node
	 */
	ForestFireGenerator(double p);

	/**
	 * Generate event stream.
	 *
	 * @param[in]	nSteps	number of time steps in the event stream
	 */
	 std::vector<GraphEvent> generate(count nSteps) override;

private:

	 double p;
	 bool firstCall;

};

} /* namespace NetworKit */

#endif /* FORESTFIREGENERATOR_H_ */

/*
 * ErdosRenyiGenerator.h
 *
 *  Created on: 21.01.2014
 *      Author: Henning
 */

#ifndef ERDOSRENYIGENERATOR_H_
#define ERDOSRENYIGENERATOR_H_

#include "StaticGraphGenerator.h"

namespace NetworKit {
/**
 * @ingroup generators
 * Creates G(n, p) graphs.
 *
 * This class is a wrapper to @ref ErdosRenyiEnumerator.
 */
class ErdosRenyiGenerator: public StaticGraphGenerator {
public:
	/**
	 * Creates random graphs in the G(n,p) model.
	 * The generation follows Vladimir Batagelj and Ulrik Brandes: "Efficient
	 * generation of large random networks", Phys Rev E 71, 036113 (2005).
	 *
	 * @param nNodes Number of nodes n in the graph.
	 * @param prob Probability of existence for each edge p.
	 * @param directed	generates a directed graph
	 */
	ErdosRenyiGenerator(count nNodes, double prob, bool directed=false);

	virtual Graph generate();

private:
	node nNodes;
	double prob;
	bool directed;
};

} /* namespace NetworKit */
#endif /* ERDOSRENYIGENERATOR_H_ */

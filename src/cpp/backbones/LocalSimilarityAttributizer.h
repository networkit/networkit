/*
 * LocalSimilarityAttributizer.h
 *
 *  Created on: 26.07.2014
 *      Author: Gerd Lindner
 */

#ifndef LOCALSIMATTRIBUTIZER_H_
#define LOCALSIMATTRIBUTIZER_H_

#include "BackboneCalculator.h"

namespace NetworKit {

/** 
 * Implementation of the Local Sparsification Algorithm by Sataluri et al.
 */
class LocalSimilarityAttributizer : public AttributeGenerator {

public:

	/**
	 * Creates a new instance of the Local Sparsification algorithm.
	 */
	LocalSimilarityAttributizer();

	EdgeAttribute getAttribute(const Graph& graph, const EdgeAttribute& attribute);

private:
	//Private helper functions
	double getSimilarity(const Graph& graph, node i, node j);
};

}
/* namespace NetworKit */
#endif /* LOCALSIMATTRIBUTIZER_H_ */

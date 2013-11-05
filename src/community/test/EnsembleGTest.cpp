/*
 * EnsembleGTest.cpp
 *
 *  Created on: 31.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#include "EnsembleGTest.h"

#include "../EnsemblePreprocessing.h"
#include "../../graph/GraphGenerator.h"
#include "../../clustering/Modularity.h"
#include "../../community/LabelPropagation.h"
#include "../../community/Louvain.h"
#include "../../overlap/HashingOverlapper.h"

namespace NetworKit {

TEST_F(EnsembleGTest, testEnsemblePreprocessing) {
	count n = 1000;
	count k = 10;
	double pin = 1.0;
	double pout = 0.0;

	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pin, pout);

	EnsemblePreprocessing ensemble;

	count b = 4;
	for (count i = 0; i < b; ++i) {
		ensemble.addBaseClusterer(*(new LabelPropagation()));
	}
	ensemble.setFinalClusterer(*(new Louvain()));
	ensemble.setOverlapper(*(new HashingOverlapper));

	Clustering zeta = ensemble.run(G);

	INFO("number of clusters:" << zeta.numberOfClusters());

	Modularity modularity;
	INFO("modularity: " << modularity.getQuality(zeta, G));



}

} /* namespace NetworKit */

#endif /*NOGTEST */

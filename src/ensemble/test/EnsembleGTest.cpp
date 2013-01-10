/*
 * EnsembleGTest.cpp
 *
 *  Created on: 31.12.2012
 *      Author: cls
 */

#include "EnsembleGTest.h"

namespace EnsembleClustering {

TEST_F(EnsembleGTest, testEnsembleClusterer) {

	EnsembleClusterer ensembleClusterer;
	// configure EnsembleClusterer
	QualityMeasure* qm = new Modularity();
	ensembleClusterer.setQualityMeasure(*qm);
	int b = 2; // number of base clusterers
	for (int i = 0; i < b; ++i) {
		Clusterer* baseClusterer = new LabelPropagation();
		ensembleClusterer.addBaseClusterer(*baseClusterer);
	}
	Clusterer* finalClusterer = new LabelPropagation();
	ensembleClusterer.setFinalClusterer(*finalClusterer);

	// generate clustered random graph with obvious community structure
	GraphGenerator graphGen;
	int64_t n = 20;
	int64_t k = 3;
	double pIn = 0.5;
	double pOut = 0.001;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pIn, pOut);


	Clustering zeta = ensembleClusterer.run(G);

	Modularity modularity;
	double mod = modularity.getQuality(zeta, G);

	INFO("modularity produced by EnsembleClusterer: " << mod);

	EXPECT_TRUE(zeta.isProper(G)) << "the resulting partition should be a proper clustering";
	EXPECT_EQ(k, zeta.numberOfClusters()) << " " << k << " clusters are easy to detect";

}

} /* namespace EnsembleClustering */

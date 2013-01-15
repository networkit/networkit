/*
 * EnsembleGTest.cpp
 *
 *  Created on: 31.12.2012
 *      Author: cls
 */

#include "EnsembleGTest.h"

namespace EnsembleClustering {

TEST_F(EnsembleGTest, testEnsembleClustererOnCliqueGraph) {

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
	// these parameters generate a clique graph
	double pIn = 1.0;
	double pOut = 0.0;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pIn, pOut);

	// DEBUG
	GraphIO graphio;
	graphio.writeAdjacencyList(G, "sandbox/G_Clique.adjlist");
	// DEBUG


	Clustering zeta = ensembleClusterer.run2(G);

	DEBUG("clustering produced by EnsembleClusterer: "); zeta.print();

	// DEBUG
	if (zeta.size() != G.numberOfNodes()) {
		ERROR("clustering produced by EnsembleClusterer has " << zeta.size() << " entries but n = " << G.numberOfNodes());
	}
	// DEBUG

	EXPECT_TRUE(zeta.isProper(G)) << "the resulting partition should be a proper clustering";
	EXPECT_EQ(k, zeta.numberOfClusters()) << " " << k << " clusters (cliques) are easy to detect";

	Modularity modularity;
	double mod = modularity.getQuality(zeta, G);
	INFO("modularity produced by EnsembleClusterer: " << mod);

}


TEST_F(EnsembleGTest, testEnsembleClustererOnAlmostCliqueGraph) {

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
	// these parameters generate a clique graph
	double pIn = 1.0;
	double pOut = 0.01;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pIn, pOut);

	// DEBUG
	GraphIO graphio;
	graphio.writeAdjacencyList(G, "sandbox/G_AlmostClique.adjlist");
	// DEBUG


	Clustering zeta = ensembleClusterer.run2(G);

	DEBUG("clustering produced by EnsembleClusterer: "); zeta.print();

	// DEBUG
	if (zeta.size() != G.numberOfNodes()) {
		ERROR("clustering produced by EnsembleClusterer has " << zeta.size() << " entries but n = " << G.numberOfNodes());
	}
	// DEBUG

	EXPECT_TRUE(zeta.isProper(G)) << "the resulting partition should be a proper clustering";
	EXPECT_EQ(k, zeta.numberOfClusters()) << " " << k << " clusters are easy to detect";

	Modularity modularity;
	double mod = modularity.getQuality(zeta, G);
	INFO("modularity produced by EnsembleClusterer: " << mod);

}

} /* namespace EnsembleClustering */

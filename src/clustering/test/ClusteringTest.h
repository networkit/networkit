/*
 * ClusteringTest.h
 *
 *  Created on: 12.12.2012
 *      Author: cls
 */

#ifndef CLUSTERINGTEST_H_
#define CLUSTERINGTEST_H_

#include <gtest/gtest.h>



#include "../Clustering.h"
#include "../Modularity.h"
#include "../ClusteringGenerator.h"
#include "../../graph/GraphGenerator.h"

namespace EnsembleClustering {

class ClusteringTest: public testing::Test {


};



TEST_F(ClusteringTest, testModularity) {
	GraphGenerator graphGenerator;

	Graph G = graphGenerator.makeCompleteGraph(20);

	ClusteringGenerator clusteringGenerator;

	Clustering singleton = clusteringGenerator.makeSingletonClustering(G);
	Clustering one = clusteringGenerator.makeOneClustering(G);

	Modularity modularity(G);

	double modSingleton = modularity.getQuality(singleton);
	double modOne = modularity.getQuality(one);

	// TODO: log values
	EXPECT_EQ(modOne, 0.0) << "1-clustering should have modularity of 0.0";
	ASSERT_LE(modSingleton, 0.0) << "singleton clustering should have modularity less than 0.0";

}

} /* namespace EnsembleClustering */
#endif /* CLUSTERINGTEST_H_ */

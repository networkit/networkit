/*
 * PropertiesGTest.cpp
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#include "PropertiesGTest.h"

namespace NetworKit {

PropertiesGTest::PropertiesGTest() {
	// TODO Auto-generated constructor stub

}

PropertiesGTest::~PropertiesGTest() {
	// TODO Auto-generated destructor stub
}


TEST_F(PropertiesGTest, testClusteringCoefficient) {

	GraphGenerator gen;
	Graph G = gen.makeErdosRenyiGraph(100, 1.0);

	ClusteringCoefficient clusteringCoefficient;
	double cc = clusteringCoefficient.calculate(G);

	EXPECT_EQ(1.0, cc);

}

TEST_F(PropertiesGTest, testDegreeDistribution) {
	Graph G(3);
	std::vector<count> degreeDist = GraphProperties::degreeDistribution(G);
	G.addEdge(0,1);
	G.addEdge(1,2);
	G.addEdge(2,0);
	degreeDist = GraphProperties::degreeDistribution(G);
	EXPECT_EQ(0, degreeDist[0]);
	EXPECT_EQ(0, degreeDist[1]);
	EXPECT_EQ(3, degreeDist[2]);


}



TEST_F(PropertiesGTest, testLocalClusteringCoefficients) {

	// Test case for a complete graph
	GraphGenerator gen;
	Graph G_complete = gen.makeCompleteGraph(4);

	std::vector<double> coefficients = GraphProperties::localClusteringCoefficients(G_complete);
	for (double cc : coefficients) {
		EXPECT_EQ(1.0, cc) << "In a clique all possible triangles are closed so all local clustering coefficients are 1";
	}

	// Test case for graph with degree-1 nodes
	Graph G1(4);
	//build a path
	G1.addEdge(0, 1);
	G1.addEdge(1, 2);
	G1.addEdge(2, 3);





}


TEST_F(PropertiesGTest, testAverageLocalClusteringCoefficient) {

	// Test case for a complete graph
	GraphGenerator gen;
	Graph G_complete = gen.makeCompleteGraph(4);

	EXPECT_EQ(1.0, GraphProperties::averageLocalClusteringCoefficient(G_complete)) << "should be 1.0 for a complete graph";

	// Test case for graph with degree-1 nodes
	Graph G_path(4);
	//build a path
	G_path.addEdge(0, 1);
	G_path.addEdge(1, 2);
	G_path.addEdge(2, 3);

	EXPECT_EQ(0.0, GraphProperties::averageLocalClusteringCoefficient(G_path)) << "should be 0.0 for a path";


}


TEST_F(PropertiesGTest, testLocalClusteringCoefficientPerDegree) {
	GraphGenerator gen;
	Graph G = gen.makeCompleteGraph(2);

	std::vector<double> coefficients = GraphProperties::localClusteringCoefficientPerDegree(G);

	EXPECT_EQ(0.0, coefficients[0]);
	EXPECT_EQ(0.0, coefficients[1]);
	EXPECT_EQ(0.0, coefficients[2]);
	EXPECT_EQ(1.0, coefficients[3]);


	// Test case for a graph, which contains nodes with degree 1

	Graph G1(5);
	G1.addEdge(0,1);

	std::vector<double> coefficients1 = GraphProperties::localClusteringCoefficientPerDegree(G1);
	for (double cc : coefficients1) {
		EXPECT_EQ(none, cc) << "Local clustering coefficients should not be calculated should not be calculated for nodes with degree 0 and 1";
	}

	G1.addEdge(1, 2);
	G1.addEdge(2, 0);

	coefficients1 = GraphProperties::localClusteringCoefficientPerDegree(G1);

	EXPECT_EQ(none, coefficients1[0]);
	EXPECT_EQ(none, coefficients1[1]);
	EXPECT_EQ(1.0, coefficients1[2]);


/*
	// Test case for a graph, which contains nodes with degree 0

	Graph G0(1);

	std::vector<double> coefficients0 = GraphProperties::localClusteringCoefficientPerDegree(G0);
	for (double cc : coefficients0) {
		EXPECT_EQ(none, cc) << "Local clustering coefficients should not be calculated should not be calculated for nodes with degree 0 and 1";
	}

	*/



}



} /* namespace NetworKit */

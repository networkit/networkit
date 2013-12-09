/*
 * PropertiesGTest.cpp
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#ifndef NOGTEST

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
	double cc = clusteringCoefficient.avgLocal(G);

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

	std::vector<double> coefficients_G1 = GraphProperties::localClusteringCoefficients(G1);

	EXPECT_EQ(0, coefficients_G1[0]);
	EXPECT_EQ(0, coefficients_G1[1]);
	EXPECT_EQ(0, coefficients_G1[2]);

	// Test case for graph with degree-1 nodes and an isolated node
	Graph G2(5);
	//build a path
	G2.addEdge(0, 1);
	G2.addEdge(1, 2);
	G2.addEdge(2, 3);
	G2.addNode();
	std::vector<double> coefficients_G2 = GraphProperties::localClusteringCoefficients(G2);

	EXPECT_EQ(0, coefficients_G2[0]);
	EXPECT_EQ(0, coefficients_G2[1]);
	EXPECT_EQ(0, coefficients_G2[2]);

}


TEST_F(PropertiesGTest, testCoreDecomposition) {
	count n = 16;
	Graph G(n);

	// create graph used in Baur et al. and network analysis lecture
	G.addEdge(2, 4);
	G.addEdge(3, 4);
	G.addEdge(4, 5);
	G.addEdge(5, 7);
	G.addEdge(6, 7);

	G.addEdge(6, 8);
	G.addEdge(6, 9);
	G.addEdge(6, 11);
	G.addEdge(7, 12);
	G.addEdge(8, 9);

	G.addEdge(8, 10);
	G.addEdge(8, 11);
	G.addEdge(8, 13);
	G.addEdge(9, 10);
	G.addEdge(9, 11);

	G.addEdge(9, 13);
	G.addEdge(10, 11);
	G.addEdge(10, 13);
	G.addEdge(10, 14);
	G.addEdge(11, 13);

	G.addEdge(11, 14);
	G.addEdge(12, 15);
	G.addEdge(13, 14);
	G.addEdge(14, 15);

	EXPECT_EQ(n, G.numberOfNodes()) << "should have " << n << " vertices";
	EXPECT_EQ(24, G.numberOfEdges()) << "should have 24 edges";

	// compute core decomposition
	CoreDecomposition coreDec;
	std::vector<count> coreness = coreDec.run(G);

	EXPECT_EQ(0, coreness[0]) << "expected coreness";
	EXPECT_EQ(0, coreness[1]) << "expected coreness";
	EXPECT_EQ(1, coreness[2]) << "expected coreness";
	EXPECT_EQ(1, coreness[3]) << "expected coreness";
	EXPECT_EQ(1, coreness[4]) << "expected coreness";
	EXPECT_EQ(1, coreness[5]) << "expected coreness";
	EXPECT_EQ(3, coreness[6]) << "expected coreness";
	EXPECT_EQ(2, coreness[7]) << "expected coreness";
	EXPECT_EQ(4, coreness[8]) << "expected coreness";
	EXPECT_EQ(4, coreness[9]) << "expected coreness";
	EXPECT_EQ(4, coreness[10]) << "expected coreness";
	EXPECT_EQ(4, coreness[11]) << "expected coreness";
	EXPECT_EQ(2, coreness[12]) << "expected coreness";
	EXPECT_EQ(4, coreness[13]) << "expected coreness";
	EXPECT_EQ(3, coreness[14]) << "expected coreness";
	EXPECT_EQ(2, coreness[15]) << "expected coreness";
}

/*
TEST_F(PropertiesGTest, testCoreDecompositionOnGraphFiles) {
	CoreDecomposition coreDec;
  METISGraphReader input;
  
  Graph G;
  std::vector<count> corenesses;
  std::ofstream output;
  
  std::list<std::string> graphList = {"celegans_metabolic", "polblogs", "hep-th"};
  
  std::for_each(graphList.begin(), graphList.end(), [&](std::string filename){
    G = input.read(std::string("input/") + filename.str() + ".graph");
    corenesses = coreDec.run(G);
    output.open(std::string("output/") + filename + ".sol");
    std::for_each(corenesses.begin(), corenesses.end(), [&](count coreness){
      output << coreness << std::endl;
    });
    output.close();
  });
}
*/


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

	// Test case for graph with degree-1 nodes
	Graph G_path_isolated(4);
	//build a path
	G_path_isolated.addEdge(0, 1);
	G_path_isolated.addEdge(1, 2);
	G_path_isolated.addEdge(2, 3);
	G_path_isolated.addNode();

	EXPECT_EQ(0.0, GraphProperties::averageLocalClusteringCoefficient(G_path_isolated)) << "should be 0.0 for a path";


}


TEST_F(PropertiesGTest, testLocalClusteringCoefficientPerDegree) {
	GraphGenerator gen;
	Graph G = gen.makeCompleteGraph(5);

	std::vector<double> coefficients = GraphProperties::localClusteringCoefficientPerDegree(G);

	EXPECT_EQ(0.0, coefficients[0]);
	EXPECT_EQ(0.0, coefficients[1]);
	EXPECT_EQ(0.0, coefficients[2]);
	EXPECT_EQ(0.0, coefficients[3]);
	EXPECT_EQ(1.0, coefficients[4]);

	// Test case for a graph, which contains nodes with degree 1

	Graph G1(5);
	G1.addEdge(0,1);

	std::vector<double> coefficients1 = GraphProperties::localClusteringCoefficientPerDegree(G1);
	for (double cc : coefficients1) {
		EXPECT_EQ(0, cc) << "Local clustering coefficients should not be calculated should not be calculated for nodes with degree 0 and 1";
	}

	G1.addEdge(1, 2);
	G1.addEdge(2, 0);
	G1.addNode();

	coefficients1 = GraphProperties::localClusteringCoefficientPerDegree(G1);

	EXPECT_EQ(0, coefficients1[0]);
	EXPECT_EQ(0, coefficients1[1]);
	EXPECT_EQ(1.0, coefficients1[2]);

}

TEST_F(PropertiesGTest, testLocalClusteringCoefficientOnARealGraph) {
	// Reading the graph
	std::string path = "input/jazz.graph";

	METISGraphReader reader;
	Graph G = reader.read(path);
	count n = 198;
	count m = 2742;
	EXPECT_FALSE(G.isEmpty());
	EXPECT_EQ(n, G.numberOfNodes()) << "There are " << n << " nodes in the  graph";
	EXPECT_EQ(m, G.numberOfEdges()) << "There are " << m << " edges in the  graph";

	// Calculating the parameters
	std::vector<count> 	degreeDist = GraphProperties::degreeDistribution(G);
	std::vector<double> coefficients = GraphProperties::localClusteringCoefficients(G);
	std::vector<double> coefficientsPerDegree = GraphProperties::localClusteringCoefficientPerDegree(G);
	double avgCoefficient = GraphProperties::averageLocalClusteringCoefficient(G);



	// Comparing the results with the values calculated by the NetworkX package (http://networkx.github.io/)

	// NetworkX: nx.degree_histogram(G)
	std::vector<count> degreeDistNetworkX = {0, 5, 3, 3, 3, 4, 5, 3, 3, 3, 2, 1, 3, 6, 4, 1, 4, 5, 4, 7, 8, 2, 2, 9, 5,
			5, 3, 5, 3, 5, 2, 9, 3, 4, 3, 1, 2, 3, 1, 6, 5, 4, 3, 3, 1, 3, 5, 0, 2, 2, 0, 1, 2, 2, 1, 2, 2, 1, 0, 2, 2,
			0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
			0, 0, 0, 1
	};


	// NetworkX: nx.clustering(G)
	std::vector<double> coefficientsNetworkX = {0.6561264822134387, 1.0, 1.0, 0.825, 0.421256038647343, 1.0, 0.4423558897243108,
			0.7428571428571429, 0.5526315789473685, 0.6780626780626781, 0.49407114624505927, 0.7285714285714285,
			0.6733333333333333, 0.7066666666666667, 0.9285714285714286, 0.8055555555555556, 0.0, 0.46546546546546547,
			0.6557971014492754, 0.5636363636363636, 0.0, 1.0, 0.7564102564102564, 0.4802955665024631, 0.6, 0.6903225806451613,
			0.9743589743589743, 0.7486772486772487, 0.6580645161290323, 0.0, 0.7210526315789474, 0.6109756097560975,
			0.7914438502673797, 0.5, 0.6356589147286822, 0.7353846153846154, 0.5, 0.5847953216374269, 0.6, 0.7914438502673797,
			1.0, 0.5151515151515151, 0.6083743842364532, 0.7447447447447447, 1.0, 0.6482213438735178, 0.6666666666666666,
			0.9285714285714286, 0.4738675958188153, 0.49002849002849, 0.7192982456140351, 0.8666666666666667, 0.5052264808362369,
			0.4473429951690821, 0.7866666666666666, 0.6491935483870968, 0.551051051051051, 0.8649193548387096, 0.6818181818181818,
			0.31162280701754386, 0.5416666666666666, 0.9333333333333333, 0.9118279569892473, 0.8891129032258065, 0.8617424242424242,
			0.8306595365418895, 0.5862068965517241, 0.8241758241758241, 0.44702467343976776, 0.39767318878900054, 0.4689655172413793,
			0.8333333333333334, 0.5612903225806452, 0.5157894736842106, 0.7720588235294118, 0.875, 0.6798029556650246, 0.46153846153846156,
			0.6210526315789474, 0.47368421052631576, 0.4831591173054588, 1.0, 0.38340151957919344, 0.5824175824175825, 0.8771929824561403,
			0.5282051282051282, 0.6666666666666666, 0.42063492063492064, 0.6838235294117647, 0.5454545454545454, 0.4746376811594203, 0.7,
			0.7663817663817664, 0.796923076923077, 0.5128205128205128, 0.3717948717948718, 0.7843137254901961, 0.6405797101449275,
			0.46497175141242936, 0.5952380952380952, 0.592687074829932, 0.8947368421052632, 0.4624505928853755, 0.5157894736842106,
			0.5718085106382979, 0.6572199730094467, 0.6699857752489331, 0.4858757062146893, 0.8352272727272727, 0.6323366555924695,
			0.43902439024390244, 0.6825396825396826, 1.0, 0.47463002114164904, 0.38095238095238093, 0.46842105263157896, 0.6794871794871795,
			0.41585365853658535, 0.7777777777777778, 0.0, 0.45209176788124156, 0.4727272727272727, 0.6333333333333333, 0.8333333333333334,
			0.6580645161290323, 0.8676470588235294, 0.6477832512315271, 0.5525641025641026, 0.6798418972332015, 0.6, 0.512987012987013,
			0.40396396396396395, 0.5333333333333333, 0.5512820512820513, 0.46901960784313723, 0.23717171717171717, 0.5833333333333334,
			0.6679841897233202, 0.6013071895424836, 0.7075268817204301, 0.7312252964426877, 0.524731182795699, 0.6176470588235294,
			0.8014705882352942, 1.0, 0.7720797720797721, 0.573549257759784, 0.7, 0.3494060097833683, 0.5222672064777328, 0.7272727272727273,
			0.6666666666666666, 0.2865853658536585, 0.6356589147286822, 0.5684210526315789, 0.8051948051948052, 0.3333333333333333,
			0.39976621858562245, 0.4666666666666667, 0.0, 0.610752688172043, 0.4152046783625731, 0.7424242424242424, 0.4702467343976778,
			0.0, 0.8241758241758241, 0.44242424242424244, 0.36208811551277303, 0.8142292490118577, 0.475177304964539, 0.4626262626262626,
			0.5294871794871795, 0.7971014492753623, 0.39064856711915535, 0.7450980392156863, 0.7318840579710145, 0.7368421052631579,
			0.5263157894736842, 0.5342995169082125, 0.4, 0.6, 0.73, 0.39541160593792174, 0.8354978354978355, 0.6666666666666666, 1.0,
			0.7647058823529411, 0.4666666666666667, 0.4924731182795699, 0.7720797720797721, 0.8066666666666666, 0.44242424242424244,
			0.7486772486772487, 0.45194805194805193, 0.4888888888888889, 0.48282828282828283, 0.6031746031746031, 0.7692307692307693
	};

	//for( std::vector<double>::const_iterator i = coefficientsNetworkX.begin(); i != coefficients.end(); ++i)
		    //std::cout << *i << ' ';

	/*for (double x : coefficients)
	{
		std::cout << x << " ";
	}*/

	EXPECT_TRUE(coefficientsNetworkX == coefficients);
	EXPECT_TRUE(degreeDistNetworkX == degreeDist);

	// Networkx: nx.average_clustering(G)
	double avgCoefficientNetworkX = 0.6174507021536301;

	EXPECT_EQ(avgCoefficientNetworkX, avgCoefficient);



}






} /* namespace NetworKit */

#endif /*NOGTEST*/

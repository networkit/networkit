/*
 * IOGTest.cpp
 *
 *  Created on: 12.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#include "IOGTest.h"


namespace NetworKit {

TEST_F(IOGTest, testGraphIOEdgeList) {
	GraphGenerator graphGen;
	Graph G = graphGen.makeCircularGraph(20);
	GraphIO graphio;
	std::string path = "output/edgelist.txt";
	graphio.writeEdgeList(G, path);

	bool exists = false;
	std::ifstream file(path);
	if (file) {
		exists = true;
	}
	EXPECT_TRUE(exists) << "A file should have been created : " << path;
}

TEST_F(IOGTest, testGraphIOAdjacencyList) {
	GraphGenerator graphGen;
	Graph G = graphGen.makeCircularGraph(20);
	GraphIO graphio;
	std::string path = "output/circular.adjlist";
	graphio.writeAdjacencyList(G, path);

	bool exists = false;
	std::ifstream file(path);
	if (file) {
		exists = true;
	}
	EXPECT_TRUE(exists) << "A file should have been created : " << path;
}


TEST_F(IOGTest, testGraphIOForIsolatedNodes) {
	GraphGenerator graphGen;
	Graph G(20);
	GraphIO graphio;
	std::string path = "output/isolated.adjlist";
	graphio.writeAdjacencyList(G, path);

	bool exists = false;
	std::ifstream file(path);
	if (file) {
		exists = true;
	}
	EXPECT_TRUE(exists) << "A file should have been created : " << path;
}



TEST_F(IOGTest, testMETISGraphReader) {
	std::string path = "input/jazz.graph";

	METISGraphReader reader;
	Graph G = reader.read(path);
	count n = 198;
	count m = 2742;
	EXPECT_FALSE(G.isEmpty());
	EXPECT_EQ(n, G.numberOfNodes()) << "There are " << n << " nodes in the  graph";
	EXPECT_EQ(m, G.numberOfEdges()) << "There are " << m << " edges in the  graph";
}

TEST_F(IOGTest, testMETISGraphReaderWithWeights) {
	std::string path = "input/lesmis.graph";

	METISGraphReader reader;
	Graph G = reader.read(path);

	EXPECT_FALSE(G.isEmpty());
	count n = 77;
	count m = 254;
	EXPECT_EQ(n, G.numberOfNodes()) << "There are " << n << " nodes in the  graph";
	EXPECT_EQ(m, G.numberOfEdges()) << "There are " << m << " edges in the  graph";
}

TEST_F(IOGTest, testMETISGraphWriter) {
	std::string path = "input/jazz1.graph";
	Graph G = Graph(3);
	G.addEdge(0,2);
	G.addEdge(1,1);
	G.addEdge(1,2);
	G.addEdge(2,2);

	METISGraphWriter writer;
	writer.write(G, false, path);
    bool exists = false;
	std::ifstream file(path);
	if (file) {
		exists = true;
	}
	EXPECT_TRUE(exists) << "A file should have been created : " << path;

}

TEST_F(IOGTest, testMETISGraphWriterWithWeights) {
	std::string path = "input/jazz2.graph";
	Graph G = Graph(5);
	G.addEdge(0,2);
	G.addEdge(0,1);
	G.addEdge(0,0);
	G.addEdge(1,1);

	METISGraphWriter writer;
	writer.write(G, true, path);
    bool exists = false;
	std::ifstream file(path);
	if (file) {
		exists = true;
	}
	EXPECT_TRUE(exists) << "A file should have been created : " << path;

}

TEST_F(IOGTest, testClusteringWriterAndReader) {
	// write clustering first
	std::string path = "output/example.clust";

	GraphGenerator graphGen;
	count n = 100;
	count k = 3;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusteringGen;
	Clustering zeta = clusteringGen.makeRandomClustering(G, k);

	ClusteringWriter writer;
	writer.write(zeta, path);

	// check if file exists
	bool exists = false;
	std::ifstream file(path);
	if (file) {
		exists = true;
	}
	EXPECT_TRUE(exists) << "clustering file should have been written to: " << path;


	ClusteringReader reader;
	Clustering read = reader.read(path);

	EXPECT_EQ(n, read.numberOfNodes()) << "read clustering should contain n nodes";
	EXPECT_TRUE(read.isProper(G)) << "read clustering should be proper clustering of G";
	EXPECT_TRUE(read.equals(zeta, G)) << "read clustering should be identical to created clustering";
}


TEST_F(IOGTest, testDotGraphWriter) {
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(42);

	std::string path = "output/example.dot";

	DotGraphWriter writer;
	writer.write(G, path);

	// check if file exists
	bool exists = false;
	std::ifstream file(path);
	if (file) {
		exists = true;
	}
	EXPECT_TRUE(exists) << "graph file should have been written to: " << path;
}

TEST_F(IOGTest, tryDGSReaderOnBigFile) {
	// read example graph
	DGSReader reader;
	Graph G;
	GraphEventProxy Gproxy(G);
	reader.read("/Users/forigem/KIT/NetworKit-CommunityDetection/input/AuthorsGraph.dgs", Gproxy);
}



TEST_F(IOGTest, tryDGSReader) {
	// read example graph
	DGSReader reader;
	Graph G;
	GraphEventProxy Gproxy(G);
	reader.read("input/example2.dgs", Gproxy);

	// get input parameters
	count nodeCount = G.numberOfNodes();
	DEBUG("Number of nodes " << nodeCount);
	EXPECT_EQ(3, nodeCount);
	count edgeCount = G.numberOfEdges();
	DEBUG("Number of edges " << edgeCount);
	EXPECT_EQ(2, edgeCount);

	G.forNodes([&](node n) {
		DEBUG("DEGREE OF NODE: " << G.degree(n) << std::endl);
	});

}

TEST_F(IOGTest, testEdgeListReader) {
	EdgeListReader reader;

	Graph G = reader.read("input/example.edgelist");

}


TEST_F(IOGTest, testMETISGraphReaderWithIsolatedNodes) {
	METISGraphReader reader;
	Graph G = reader.read("input/example.graph");
	EXPECT_EQ(4, G.numberOfNodes());
	EXPECT_EQ(1, G.numberOfEdges());
	EXPECT_TRUE(G.hasNode(0));
	EXPECT_TRUE(G.hasNode(1));
	EXPECT_TRUE(G.hasNode(2));
	EXPECT_TRUE(G.hasNode(3));
}

TEST_F(IOGTest, testMETISGraphReaderForNodeExistence2) {
	METISGraphReader reader;
	Graph G = reader.read("input/jazz.graph");
	EXPECT_TRUE(G.hasNode(0));
	EXPECT_EQ(198, G.numberOfNodes());
	EXPECT_EQ(2742, G.numberOfEdges());
}



} /* namespace NetworKit */

#endif /* NOGTEST */

/*
 * APSPGTest.cpp
 *
 *  Created on: 07.07.2015
 *      Author: Arie Slobbe, Elisabetta Bergamini
 */

#ifndef NOGTEST

#include "APSPGTest.h"
#include "../APSP.h"
#include "../DynAPSP.h"
#include "../../dynamics/GraphEvent.h"
#include <string>


namespace NetworKit {

TEST_F(APSPGTest, testAPSP) {
/* Graph:
     ______
		/      \
	 0    3   6
		\  / \ /
		 2    5
		/  \ / \
	 1    4   7
*/
	int n = 8;
	Graph G(n);

	G.addEdge(0, 2);
	G.addEdge(1, 2);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 5);
	G.addEdge(4, 5);
	G.addEdge(5, 6);
	G.addEdge(5, 7);
	G.addEdge(0, 6);

	APSP apsp(G);
	apsp.run();
	std::vector<std::vector<edgeweight> > distances = apsp.getDistances();
	//INFO("distances[0]: ", distances[0][0], distances[0][1], distances[0][2], distances[0][3], distances[0][4], distances[0][5], distances[0][6]);
	//INFO("distances[1]: ", distances[1][0], distances[1][1], distances[1][2], distances[1][3], distances[1][4], distances[1][5], distances[1][6]);
}

TEST_F(APSPGTest, testAPSPInsertionUndirectedUnweighted) {
	// build G'
	Graph G(7);
	G.addEdge(0,1);
	G.addEdge(1,2);
	G.addEdge(1,3);
	G.addEdge(3,4);
	G.addEdge(4,5);
	G.addEdge(5,6);

	// run dyn apsp with insertion of edge (2, 5)
	DynAPSP apsp(G);
	apsp.run();
	std::vector<std::vector<edgeweight> > distances = apsp.getDistances();
	//INFO("distances[0]: ", distances[0][0], distances[0][1], distances[0][2], distances[0][3], distances[0][4], distances[0][5], distances[0][6]);

	// apply graph update
	G.addEdge(2, 5);
	GraphEvent event(GraphEvent::EDGE_ADDITION, 2, 5);
	apsp.update(event);
	distances = apsp.getDistances();

	apsp.run();
	std::vector<std::vector<edgeweight> > distances2 = apsp.getDistances();
	G.forNodes([&](node i) {
		G.forNodes([&](node j) {
			//INFO("i, j = ", i, j);
			EXPECT_EQ(distances[i][j], distances2[i][j]);
		});
	});

	// EXPECT_EQ(distances[0][5], 3);
	// EXPECT_EQ(distances[0][6], 4);
}

TEST_F(APSPGTest, testAPSPDeletionUndirectedUnweighted) {
	// build G'
	Graph G(7);
	G.addEdge(0,1);
	G.addEdge(1,2);
	G.addEdge(1,3);
	G.addEdge(2,5);
	G.addEdge(3,4);
	G.addEdge(4,5);
	G.addEdge(5,6);


	// run dyn apsp with deletion of edge (2,5)
	DynAPSP apsp(G);
	apsp.run();
	std::vector<std::vector<edgeweight> > distances = apsp.getDistances();
	//INFO("distances[0]: ", distances[0][0], distances[0][1], distances[0][2], distances[0][3], distances[0][4], distances[0][5], distances[0][6]);

	// apply graph update
	G.removeEdge(2, 5);
	GraphEvent event(GraphEvent::EDGE_REMOVAL, 2, 5);
	apsp.update(event);
	distances = apsp.getDistances();

	apsp.run();
	std::vector<std::vector<edgeweight> > distances2 = apsp.getDistances();
	G.forNodes([&](node i) {
		G.forNodes([&](node j) {
			//INFO("i, j = ", i, j);
			EXPECT_EQ(distances[i][j], distances2[i][j]);
		});
	});

	// EXPECT_EQ(distances[0][5], 4);
	// EXPECT_EQ(distances[0][6], 5);
}

TEST_F(APSPGTest, testAPSPInsertionDirectedUnweighted) {
	Graph G(6, false, true);
	G.addEdge(0,1);
	G.addEdge(1,2);
	G.addEdge(2,3);
	G.addEdge(3,4);
	G.addEdge(4,5);
	G.addEdge(5,0);
	DynAPSP apsp(G);
	apsp.run();
	std::vector<std::vector<edgeweight> > distances = apsp.getDistances();

	G.forNodes([&](node i) {
		std::string Output = " ";
		G.forNodes([&](node j) {
			Output += "distance[" + std::to_string(i) + "][" + std::to_string(j) + "] = " + std::to_string(distances[i][j]);
		});
		Output += "/n";
		INFO(Output);
	});

	G.addEdge(0,3);
	GraphEvent event(GraphEvent::EDGE_ADDITION, 0, 3);
	apsp.update(event);
	distances = apsp.getDistances();

	DynAPSP apsp2(G);
	apsp2.run();
	std::vector<std::vector<edgeweight> > distances2 = apsp2.getDistances();

	G.forNodes([&](node i) {
		G.forNodes([&](node j) {
			//INFO("i, j: ", i, j);
			EXPECT_EQ(distances[i][j], distances2[i][j]);
		});
	});
}

TEST_F(APSPGTest, testAPSPDeletionDirectedUnweighted) {
	Graph G(6, false, true);
	G.addEdge(0,1);
	G.addEdge(0,3);
	G.addEdge(1,2);
	G.addEdge(2,3);
	G.addEdge(3,4);
	G.addEdge(4,5);
	G.addEdge(5,0);

	DynAPSP apsp(G);
	apsp.run();
	std::vector<std::vector<edgeweight> > distances = apsp.getDistances();

	G.forNodes([&](node i) {
		std::string Output = " ";
		G.forNodes([&](node j) {
			Output += "distance[" + std::to_string(i) + "][" + std::to_string(j) + "] = " + std::to_string(distances[i][j]);
		});
		INFO(Output);
	});

	INFO("START UPDATE");

	G.removeEdge(0,3);
	GraphEvent event(GraphEvent::EDGE_REMOVAL, 0, 3);
	apsp.update(event);
	distances = apsp.getDistances();

	G.forNodes([&](node i) {
		std::string Output = " ";
		G.forNodes([&](node j) {
			Output += "distance[" + std::to_string(i) + "][" + std::to_string(j) + "] = " + std::to_string(distances[i][j]);
		});
		INFO(Output);
	});

	DynAPSP apsp2(G);
	apsp2.run();
	std::vector<std::vector<edgeweight> > distances2 = apsp2.getDistances();

	G.forNodes([&](node i) {
		G.forNodes([&](node j) {
			//INFO("i, j: ", i, j);
			EXPECT_EQ(distances[i][j], distances2[i][j]);
		});
	});
}

TEST_F(APSPGTest, testAPSPInsertionUndirectedWeighted) {
	// here
}

TEST_F(APSPGTest, testAPSPDeletionUndirectedWeighted) {
	// here
}

} /* namespace NetworKit */

#endif /*NOGTEST */

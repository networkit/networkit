/*
 * APSPGTest.cpp
 *
 *  Created on: 07.07.2015
 *      Author: Arie Slobbe
 */

#ifndef NOGTEST

#include "APSPGTest.h"
#include "../APSP.h"
#include "../DynAPSP.h"
#include "../../dynamics/GraphEvent.h"


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
	INFO("distances[0]: ", distances[0][0], distances[0][1], distances[0][2], distances[0][3], distances[0][4], distances[0][5], distances[0][6]);
	INFO("distances[1]: ", distances[1][0], distances[1][1], distances[1][2], distances[1][3], distances[1][4], distances[1][5], distances[1][6]);
}

TEST_F(APSPGTest, testTempAPSPInsertion) {
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
	INFO("distances[0]: ", distances[0][0], distances[0][1], distances[0][2], distances[0][3], distances[0][4], distances[0][5], distances[0][6]);

	// apply graph update
	G.addEdge(2, 5);
	GraphEvent event;
	event.type = GraphEvent::EDGE_ADDITION;
	event.u = 2;
	event.v = 5;
	apsp.update(event);
	distances = apsp.getDistances();
	INFO("distances[0][0] after update: ", distances[0][0]);
	// EXPECT_EQ(D[0][5], 3);
	// EXPECT_EQ(D[0][6], 4);
}

TEST_F(APSPGTest, testTempAPSPDeletion) {
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
	INFO("distances[0]: ", distances[0][0], distances[0][1], distances[0][2], distances[0][3], distances[0][4], distances[0][5], distances[0][6]);

	// apply graph update
	G.removeEdge(2, 5);
	GraphEvent event;
	event.type = GraphEvent::EDGE_REMOVAL;
	event.u = 2;
	event.v = 5;
	apsp.update(event);
	distances = apsp.getDistances();
	INFO("distances[0][0] after update: ", distances[0][0]);

	// EXPECT_EQ(D[0][5], 4);
	// EXPECT_EQ(D[0][6], 5);
}

} /* namespace NetworKit */

#endif /*NOGTEST */

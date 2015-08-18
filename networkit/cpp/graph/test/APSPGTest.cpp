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
	// build G
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
	std::vector<GraphEvent> batch;
	batch.push_back(event);
	apsp.update(batch);
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
	// build G
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
	std::vector<GraphEvent> batch;
	batch.push_back(event);
	apsp.update(batch);
	distances = apsp.getDistances();

	apsp.run();
	std::vector<std::vector<edgeweight> > distances2 = apsp.getDistances();
	G.forNodes([&](node i) {
		G.forNodes([&](node j) {
			INFO("i, j = ", i, j);
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
	std::vector<GraphEvent> batch;
	batch.push_back(event);
	apsp.update(batch);

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

	G.removeEdge(0,3);
	GraphEvent event(GraphEvent::EDGE_REMOVAL, 0, 3);
	std::vector<GraphEvent> batch;
	batch.push_back(event);
	apsp.update(batch);

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

TEST_F(APSPGTest, testAPSPUndirectedWeighted) {
	// build G
	Graph G(7, true, false);
	G.addEdge(0,1, 1);
	G.addEdge(1,2, 0.01);
	G.addEdge(1,3, 0.1);
	G.addEdge(3,4, 0.001);
	G.addEdge(4,5, 0.0001);
	G.addEdge(5,6, 0.00001);

	// Run baseline apsp with ID 1
	DynAPSP apsp(G);
	apsp.run();
	std::vector<std::vector<edgeweight> > distances = apsp.getDistances();
	INFO("distances[0]: ", distances[0][0], distances[0][1], distances[0][2], distances[0][3], distances[0][4], distances[0][5], distances[0][6]);

	// apply graph update edge addition with ID 2
	G.addEdge(2, 5, 0.000002);
	GraphEvent event2(GraphEvent::EDGE_ADDITION, 2, 5, 0.000002);
	std::vector<GraphEvent> batch2;
	batch2.push_back(event2);
	apsp.update(batch2);

	distances = apsp.getDistances();

	DynAPSP apsp2(G);
	apsp2.run();
	std::vector<std::vector<edgeweight> > distances2 = apsp2.getDistances();
	G.forNodes([&](node i) {
		G.forNodes([&](node j) {
			//INFO("i, j = ", i, j);
			EXPECT_NEAR(distances[i][j], distances2[i][j], 0.0001);
		});
	});

	// apply graph update edge weight increment with ID 3
	G.setWeight(2, 5, 0.000003);
	GraphEvent event3(GraphEvent::EDGE_WEIGHT_INCREMENT, 2, 5, 0.000001);
	std::vector<GraphEvent> batch3;
	batch3.push_back(event3);
	apsp.update(batch3);
	distances = apsp.getDistances();

	DynAPSP apsp3(G);
	apsp3.run();
	std::vector<std::vector<edgeweight> > distances3 = apsp3.getDistances();
	G.forNodes([&](node i) {
		G.forNodes([&](node j) {
			//INFO("i, j = ", i, j);
			EXPECT_NEAR(distances[i][j], distances3[i][j], 0.0001);
		});
	});

	// apply graph update edge weight update with ID 4
	G.setWeight(2, 5, 0.000001);
	GraphEvent event4(GraphEvent::EDGE_WEIGHT_UPDATE, 2, 5, 0.000001);
	std::vector<GraphEvent> batch4;
	batch4.push_back(event4);
	apsp.update(batch4);

	distances = apsp.getDistances();

	DynAPSP apsp4(G);
	apsp4.run();
	std::vector<std::vector<edgeweight> > distances4 = apsp4.getDistances();
	G.forNodes([&](node i) {
		G.forNodes([&](node j) {
			//INFO("i, j = ", i, j);
			EXPECT_NEAR(distances[i][j], distances4[i][j], 0.0001);
		});
	});

	// build G
	Graph H(7, true, false);
	H.addEdge(0,1, 1.1);
	H.addEdge(1,2, 0.1);
	H.addEdge(1,3, 0.1);
	H.addEdge(2,5, 0.1);
	H.addEdge(3,4, 0.1);
	H.addEdge(4,5, 0.1);
	H.addEdge(5,6, 0.1);

	// Run baseline apsp with ID 1
	DynAPSP Hapsp(H);
	Hapsp.run();
	std::vector<std::vector<edgeweight> > Hdistances = Hapsp.getDistances();
	INFO("distances[0]: ", distances[0][0], " ", distances[0][1]," ", distances[0][2]," ", distances[0][3]," ", distances[0][4], " ",distances[0][5]," ", distances[0][6]);

	// apply graph update edge deletion update with ID 5
	INFO("entering update 5");
	H.removeEdge(2, 5);
	GraphEvent event5(GraphEvent::EDGE_REMOVAL, 2, 5);
	std::vector<GraphEvent> batch5;
	batch5.push_back(event5);
	apsp.update(batch5);

	Hdistances = Hapsp.getDistances();

	DynAPSP apsp5(G);
	apsp5.run();
	std::vector<std::vector<edgeweight> > distances5 = apsp5.getDistances();
	G.forNodes([&](node i) {
		G.forNodes([&](node j) {
			INFO("i, j = ", i, " ", j);
			EXPECT_NEAR(distances[i][j], distances5[i][j], 0.0001);
		});
	});


}

TEST_F(APSPGTest, testAPSPUndirectedWeightedDisconnectComponent) {
	// build G
	Graph G(7, true, true);
	G.addEdge(0,1, 1.1);
	G.addEdge(1,2, 0.1);
	G.addEdge(1,3, 0.1);
	G.addEdge(3,4, 0.1);
	G.addEdge(4,5, 0.1);
	G.addEdge(5,6, 0.1);

	// Run baseline apsp with ID 1
	DynAPSP apsp(G);
	apsp.run();
	std::vector<std::vector<edgeweight> > distances = apsp.getDistances();
	INFO("distances[0]: ", distances[0][0], " ", distances[0][1]," ", distances[0][2]," ", distances[0][3]," ", distances[0][4], " ",distances[0][5]," ", distances[0][6]);

	// apply graph update edge deletion update with ID 6
	INFO("entering update 5");
	G.removeEdge(4, 5);
	GraphEvent event5(GraphEvent::EDGE_REMOVAL, 4, 5);
	std::vector<GraphEvent> batch5;
	batch5.push_back(event5);
	apsp.update(batch5);

	distances = apsp.getDistances();

	DynAPSP apsp5(G);
	apsp5.run();
	std::vector<std::vector<edgeweight> > distances5 = apsp5.getDistances();
	G.forNodes([&](node i) {
		G.forNodes([&](node j) {
			INFO("i, j = ", i, " ", j);
			EXPECT_NEAR(distances[i][j], distances5[i][j], 0.0001);
		});
	});
}

TEST_F(APSPGTest, testAPSPDirectedWeighted) {
	Graph G(5, true, true); // G+ Ghouse
	G.addEdge(3,1,1);
	G.addEdge(1,0,2);
	G.addEdge(0,2,3);
	G.addEdge(2,1,4);
	G.addEdge(1,4,5);
	G.addEdge(4,3,6);
	G.addEdge(3,2,7);
	G.addEdge(2,4,8);

	// Run baseline apsp with ID 0
	DynAPSP apsp(G);
	apsp.run();
	std::vector<std::vector<edgeweight> > distances = apsp.getDistances();
	INFO("distances[3]: ", distances[3][0], " ", distances[3][1]," ", distances[3][2]," ", distances[3][3]," ", distances[3][4]);

	// apply graph update edge deletion update with ID 1
	INFO("entering update 1");
	G.removeEdge(3, 1);
	GraphEvent event(GraphEvent::EDGE_REMOVAL, 3, 1);
	std::vector<GraphEvent> batch;
	batch.push_back(event);
	apsp.update(batch);

	distances = apsp.getDistances();

	DynAPSP apsp1(G);
	apsp1.run();
	std::vector<std::vector<edgeweight> > distances1 = apsp1.getDistances();
	G.forNodes([&](node i) {
		G.forNodes([&](node j) {
			INFO("i, j = ", i, " ", j);
			EXPECT_NEAR(distances[i][j], distances1[i][j], 0.0001);
		});
	});

	// apply graph update edge insertion update with ID 2
	INFO("entering update 2");
	G.addEdge(3, 1, 1);
	GraphEvent event2(GraphEvent::EDGE_ADDITION, 3, 1, 1);
	std::vector<GraphEvent> batch2;
	batch2.push_back(event2);
	apsp.update(batch2);

	distances = apsp.getDistances();

	DynAPSP apsp2(G);
	apsp2.run();
	std::vector<std::vector<edgeweight> > distances2 = apsp2.getDistances();
	G.forNodes([&](node i) {
		G.forNodes([&](node j) {
			INFO("i, j = ", i, " ", j);
			EXPECT_NEAR(distances[i][j], distances2[i][j], 0.0001);
		});
	});
}

TEST_F(APSPGTest, testAPSPBatchDelIns) {
	Graph G(4, true, true);
	G.addEdge(0, 1, 2);
	G.addEdge(0, 3, 2);
	G.addEdge(1, 2, 3);
	G.addEdge(2, 3, 1);

	// Run baseline apsp with ID 0
	DynAPSP apsp(G);
	apsp.run();
	std::vector<std::vector<edgeweight> > distances = apsp.getDistances();
	INFO("distances[0]: ", distances[0][0], " ", distances[0][1]," ", distances[0][2]," ", distances[0][3]);

	// apply batch update
	G.removeEdge(0, 3);
	G.addEdge(0,2, 1);
	GraphEvent event1(GraphEvent::EDGE_REMOVAL, 0, 3);
	GraphEvent event2(GraphEvent::EDGE_ADDITION, 0, 2, 1);
	std::vector<GraphEvent> batch;
	batch.push_back(event1);
	batch.push_back(event2);
	apsp.update(batch);
	distances = apsp.getDistances();

	DynAPSP apsp1(G);
	apsp1.run();
	std::vector<std::vector<edgeweight> > distances1 = apsp1.getDistances();
	G.forNodes([&](node i) {
		G.forNodes([&](node j) {
			INFO("i, j = ", i, " ", j);
			EXPECT_NEAR(distances[i][j], distances1[i][j], 0.0001);
		});
	});


}

TEST_F(APSPGTest, testAPSPBatchInsDel) {
	Graph G(5, true, true);
	G.addEdge(0, 1, 1);
	G.addEdge(0, 3, 1);
	G.addEdge(1, 2, 1);
	G.addEdge(2, 4, 1);

	// Run baseline apsp with ID 0
	DynAPSP apsp(G);
	apsp.run();
	std::vector<std::vector<edgeweight> > distances = apsp.getDistances();
	INFO("distances[0]: ", distances[0][0], " ", distances[0][1]," ", distances[0][2]," ", distances[0][3], " ", distances[0][4]);

	// apply batch update
	G.addEdge(3, 4, 1);
	G.removeEdge(0, 3);
	GraphEvent event1(GraphEvent::EDGE_ADDITION, 3, 4, 1);
	GraphEvent event2(GraphEvent::EDGE_REMOVAL, 0, 3);
	std::vector<GraphEvent> batch;
	batch.push_back(event1);
	batch.push_back(event2);
	apsp.update(batch);
	distances = apsp.getDistances();

	DynAPSP apsp1(G);
	apsp1.run();
	std::vector<std::vector<edgeweight> > distances1 = apsp1.getDistances();
	G.forNodes([&](node i) {
		G.forNodes([&](node j) {
			INFO("i, j = ", i, " ", j);
			EXPECT_EQ(distances[i][j], distances1[i][j]);
		});
	});


}

} /* namespace NetworKit */

#endif /*NOGTEST */

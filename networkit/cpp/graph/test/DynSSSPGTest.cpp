/*
 * dynSSSPGTest.cpp
 *
 *  Created on: 21.07.2014
 *      Author: ebergamini
 */

#include "DynSSSPGTest.h"
#include "../DynBFS.h"
#include "../BFS.h"
#include "../DynDijkstra.h"
#include "../Dijkstra.h"
#include "../../io/METISGraphReader.h"
#include "../../auxiliary/Log.h"


namespace NetworKit {

TEST_F(DynSSSPGTest, testDynamicBFS_1edge) {
/* Graph:
    0    3   6
     \  / \ /
      2    5
     /  \ / \
    1    4   7
 */
	count n = 8;
	Graph G(n);

	G.addEdge(0, 2);
	G.addEdge(1, 2);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 5);
	G.addEdge(4, 5);
	G.addEdge(5, 6);
	G.addEdge(5, 7);

  BFS bfs(G, 0);
  bfs.run();
  DynBFS dbfs(G, 0);
	dbfs.init();
	std::vector<GraphEvent> batch(1);
	batch[0].type = GraphEvent::EDGE_ADDITION;
	batch[0].u = 0;
	batch[0].v = 6;
	batch[0].w = 1.0;
	for (GraphEvent edge : batch) {
		G.addEdge(edge.u, edge.v, edge.w);
	}
	dbfs.update(batch);
	bfs.run();
	G.forNodes([&] (node i) {
		EXPECT_EQ(bfs.distance(i), dbfs.distance(i));
		EXPECT_EQ(bfs.numberOfPaths(i), dbfs.numberOfPaths(i));
	});
}

TEST_F(DynSSSPGTest, testDynamicBFS_batch) {
/* Graph:
		0    3   6
		\  / \ /
			2    5
		/  \ / \
		1    4   7
*/
	count n = 8;
	Graph G(n);

	G.addEdge(0, 2);
	G.addEdge(1, 2);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 5);
	G.addEdge(4, 5);
	G.addEdge(5, 6);
	G.addEdge(5, 7);

	BFS bfs(G, 0);
	bfs.run();
	DynBFS dbfs(G, 0);
	dbfs.init();
	std::vector<GraphEvent> batch(3);
	batch[0].type = GraphEvent::EDGE_ADDITION;
	batch[0].u = 3;
	batch[0].v = 7;
	batch[0].w = 1.0;
	batch[1].type = GraphEvent::EDGE_ADDITION;
	batch[1].u = 0;
	batch[1].v = 5;
	batch[1].w = 1.0;
	batch[2].type = GraphEvent::EDGE_ADDITION;
	batch[2].u = 2;
	batch[2].v = 7;
	batch[2].w = 1.0;
	for (GraphEvent edge : batch) {
		G.addEdge(edge.u, edge.v, edge.w);
	}
	dbfs.update(batch);
	bfs.run();
	G.forNodes([&] (node i) {
		EXPECT_EQ(bfs.distance(i), dbfs.distance(i));
		EXPECT_EQ(bfs.numberOfPaths(i), dbfs.numberOfPaths(i));
	});

}


TEST_F(DynSSSPGTest, testDynamicDijkstra) {
 /* Graph:
    0    3   6
     \  / \ /
      2 -- 5
     /  \ / \
    1    4   7

    Edges in the upper row have weight 3,
    the edge in the middle row has weight 1.5,
    edges in the lower row have weight 2.
 */
	count n = 8;
	Graph G(n, true);

	G.addEdge(0, 2, 3);
	G.addEdge(1, 2, 2);
	G.addEdge(2, 3, 3);
	G.addEdge(2, 4, 2);
	G.addEdge(2, 5, 1.5);
	G.addEdge(3, 5, 3);
	G.addEdge(4, 5, 2);
	G.addEdge(5, 6, 3);
	G.addEdge(5, 7, 2);

	Dijkstra dij(G, 0);
	dij.run();
	DynDijkstra ddij(G, 0);
	ddij.run();
	std::vector<GraphEvent> batch(3);
	batch[0].type = GraphEvent::EDGE_ADDITION;
	batch[0].u = 0;
	batch[0].v = 4;
	batch[0].w = 1.0;
	batch[1].type = GraphEvent::EDGE_ADDITION;
	batch[1].u = 1;
	batch[1].v = 4;
	batch[1].w = 1.0;
	batch[2].type = GraphEvent::EDGE_ADDITION;
	batch[2].u = 6;
	batch[2].v = 7;
	batch[2].w = 3.0;
	for (GraphEvent edge : batch) {
		G.addEdge(edge.u, edge.v, edge.w);
	}
	ddij.update(batch);
	dij.run();
	G.forNodes([&] (node i) {
		EXPECT_EQ(dij.distance(i), ddij.distance(i));
		EXPECT_EQ(dij.numberOfPaths(i), ddij.numberOfPaths(i));
	});

}

} /* namespace NetworKit */

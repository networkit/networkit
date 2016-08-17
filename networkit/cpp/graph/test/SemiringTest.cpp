/*
 * SemiringTest.cpp
 *
 *  Created on: 21.07.2014
 *      Author: ebergamini
 */

#include "SemiringTest.h"
#include "../DynBFS.h"
#include "../BFS.h"
#include "../DynDijkstra.h"
#include "../GenericDijkstra.h"
#include "../GenericFloydWarshall.h"
#include "../../io/METISGraphReader.h"
#include "../../auxiliary/Log.h"

#include <stack>

#include "../ShortestPathSemiring.h"
#include "../SafestPathSemiring.h"
#include "../TimeVaryingSemiring.h"


namespace NetworKit {

TEST_F(SemiringTest, testDijkstraShortestPath) {

	int n = 6;
	Graph G(n, true);

	G.addEdge(0, 1);
	G.addEdge(0, 2);
	G.addEdge(1, 3);
	G.addEdge(2, 3);
	G.addEdge(3, 4);
	G.addEdge(3, 5);
    G.indexEdges();

    std::vector<ShortestPathSemiring> SREdges(6);
    SREdges[0] = 1;
    SREdges[1] = 1;
    SREdges[2] = 1;
    SREdges[3] = 2;
    SREdges[4] = 1;
    SREdges[5] = 2;

	GenericDijkstra<ShortestPathSemiring> sssp(G, 0, SREdges);

    sssp.run();
    EXPECT_EQ(sssp.distance(0), 0);
    EXPECT_EQ(sssp.distance(1), 1);
    EXPECT_EQ(sssp.distance(2), 1);
    EXPECT_EQ(sssp.distance(3), 2);
    EXPECT_EQ(sssp.distance(4), 3);
    EXPECT_EQ(sssp.distance(5), 4);
    
};

TEST_F(SemiringTest, testDijkstraSafestPath) {

	int n = 6;
	Graph G(n, true);

	G.addEdge(0, 1);
	G.addEdge(0, 2);
	G.addEdge(1, 3);
	G.addEdge(2, 3);
	G.addEdge(3, 4);
	G.addEdge(3, 5);
    G.indexEdges();

    std::vector<SafestPathSemiring> SREdges(6);
    SREdges[0] = 0.5;
    SREdges[1] = 0.5;
    SREdges[2] = 0.7;
    SREdges[3] = 0.5;
    SREdges[4] = 0.3;
    SREdges[5] = 0.6;

	GenericDijkstra<SafestPathSemiring> sssp(G, 0, SREdges);

    sssp.run();
    EXPECT_EQ(sssp.distance(0), 1);
    EXPECT_EQ(sssp.distance(1), 0.5);
    EXPECT_EQ(sssp.distance(2), 0.5);
    EXPECT_EQ(sssp.distance(3), 0.35);
    EXPECT_EQ(sssp.distance(4), 0.105);
    EXPECT_EQ(sssp.distance(5), 0.21);
};

TEST_F(SemiringTest, testFloydWarshallShortestPath) {

	int n = 6;
	Graph G(n, true);

	G.addEdge(0, 1);
	G.addEdge(0, 2);
	G.addEdge(1, 3);
	G.addEdge(2, 3);
	G.addEdge(3, 4);
	G.addEdge(3, 5);
    G.indexEdges();

    std::vector<ShortestPathSemiring> SREdges(6);
    SREdges[0] = 1;
    SREdges[1] = 1;
    SREdges[2] = 1;
    SREdges[3] = 2;
    SREdges[4] = 1;
    SREdges[5] = 2;

	GenericFloydWarshall<ShortestPathSemiring> apsp(G, SREdges);

    apsp.run();
    EXPECT_EQ(apsp.getDistance(0, 0), 0);
    EXPECT_EQ(apsp.getDistance(0, 1), 1);
    EXPECT_EQ(apsp.getDistance(0, 2), 1);
    EXPECT_EQ(apsp.getDistance(0, 3), 2);
    EXPECT_EQ(apsp.getDistance(0, 4), 3);
    EXPECT_EQ(apsp.getDistance(0, 5), 4);
    EXPECT_EQ(apsp.getDistance(1, 0), 1);
    EXPECT_EQ(apsp.getDistance(1, 1), 0);
    EXPECT_EQ(apsp.getDistance(1, 2), 2);
    EXPECT_EQ(apsp.getDistance(1, 3), 1);
    EXPECT_EQ(apsp.getDistance(1, 4), 2);
    EXPECT_EQ(apsp.getDistance(1, 5), 3);
    EXPECT_EQ(apsp.getDistance(2, 0), 1);
    EXPECT_EQ(apsp.getDistance(2, 1), 2);
    EXPECT_EQ(apsp.getDistance(2, 2), 0);
    EXPECT_EQ(apsp.getDistance(2, 3), 2);
    EXPECT_EQ(apsp.getDistance(2, 4), 3);
    EXPECT_EQ(apsp.getDistance(2, 5), 4);
    EXPECT_EQ(apsp.getDistance(3, 0), 2);
    EXPECT_EQ(apsp.getDistance(3, 1), 1);
    EXPECT_EQ(apsp.getDistance(3, 2), 2);
    EXPECT_EQ(apsp.getDistance(3, 3), 0);
    EXPECT_EQ(apsp.getDistance(3, 4), 1);
    EXPECT_EQ(apsp.getDistance(3, 5), 2);
    EXPECT_EQ(apsp.getDistance(4, 0), 3);
    EXPECT_EQ(apsp.getDistance(4, 1), 2);
    EXPECT_EQ(apsp.getDistance(4, 2), 3);
    EXPECT_EQ(apsp.getDistance(4, 3), 1);
    EXPECT_EQ(apsp.getDistance(4, 4), 0);
    EXPECT_EQ(apsp.getDistance(4, 5), 3);
    EXPECT_EQ(apsp.getDistance(5, 0), 4);
    EXPECT_EQ(apsp.getDistance(5, 1), 3);
    EXPECT_EQ(apsp.getDistance(5, 2), 4);
    EXPECT_EQ(apsp.getDistance(5, 3), 2);
    EXPECT_EQ(apsp.getDistance(5, 4), 3);
    EXPECT_EQ(apsp.getDistance(5, 5), 0);
};

TEST_F(SemiringTest, testFloydWarshallSafestPath) {

	int n = 6;
	Graph G(n, true);

	G.addEdge(0, 1);
	G.addEdge(0, 2);
	G.addEdge(1, 3);
	G.addEdge(2, 3);
	G.addEdge(3, 4);
	G.addEdge(3, 5);
    G.indexEdges();

    std::vector<SafestPathSemiring> SREdges(6);
    SREdges[0] = 0.5;
    SREdges[1] = 0.5;
    SREdges[2] = 0.7;
    SREdges[3] = 0.5;
    SREdges[4] = 0.3;
    SREdges[5] = 0.6;

	GenericFloydWarshall<SafestPathSemiring> apsp(G, SREdges);

    apsp.run();
    EXPECT_EQ(apsp.getDistance(0, 0), 1);
    EXPECT_EQ(apsp.getDistance(0, 1), 0.5);
    EXPECT_EQ(apsp.getDistance(0, 2), 0.5);
    EXPECT_EQ(apsp.getDistance(0, 3), 0.35);
    EXPECT_EQ(apsp.getDistance(0, 4), 0.105);
    EXPECT_EQ(apsp.getDistance(0, 5), 0.21);
    EXPECT_EQ(apsp.getDistance(1, 0), 0.5);
    EXPECT_EQ(apsp.getDistance(1, 1), 1);
    EXPECT_EQ(apsp.getDistance(1, 2), 0.35);
    EXPECT_EQ(apsp.getDistance(1, 3), 0.7);
    EXPECT_EQ(apsp.getDistance(1, 4), 0.21);
    EXPECT_EQ(apsp.getDistance(1, 5), 0.42);
    EXPECT_EQ(apsp.getDistance(2, 0), 0.5);
    EXPECT_EQ(apsp.getDistance(2, 1), 0.35);
    EXPECT_EQ(apsp.getDistance(2, 2), 1);
    EXPECT_EQ(apsp.getDistance(2, 3), 0.5);
    EXPECT_EQ(apsp.getDistance(2, 4), 0.15);
    EXPECT_EQ(apsp.getDistance(2, 5), 0.3);
    EXPECT_EQ(apsp.getDistance(3, 0), 0.35);
    EXPECT_EQ(apsp.getDistance(3, 1), 0.7);
    EXPECT_EQ(apsp.getDistance(3, 2), 0.5);
    EXPECT_EQ(apsp.getDistance(3, 3), 1);
    EXPECT_EQ(apsp.getDistance(3, 4), 0.3);
    EXPECT_EQ(apsp.getDistance(3, 5), 0.6);
    EXPECT_EQ(apsp.getDistance(4, 0), 0.105);
    EXPECT_EQ(apsp.getDistance(4, 1), 0.21);
    EXPECT_EQ(apsp.getDistance(4, 2), 0.15);
    EXPECT_EQ(apsp.getDistance(4, 3), 0.3);
    EXPECT_EQ(apsp.getDistance(4, 4), 1);
    EXPECT_EQ(apsp.getDistance(4, 5), 0.18);
    EXPECT_EQ(apsp.getDistance(5, 0), 0.21);
    EXPECT_EQ(apsp.getDistance(5, 1), 0.42);
    EXPECT_EQ(apsp.getDistance(5, 2), 0.3);
    EXPECT_EQ(apsp.getDistance(5, 3), 0.6);
    EXPECT_EQ(apsp.getDistance(5, 4), 0.18);
    EXPECT_EQ(apsp.getDistance(5, 5), 1);
};

TEST_F(SemiringTest, testFloydWarshallTimeVarying) {

	int n = 5;
	Graph G(n, true, true);

	G.addEdge(0, 1);
	G.addEdge(1, 2);
	G.addEdge(1, 3);
	G.addEdge(2, 3);
	G.addEdge(3, 4);
    G.indexEdges();

    std::vector<TimeVaryingSemiring> SREdges(5);
    SREdges[0] = TimeVaryingSemiring({{3,3,4,4}});
    SREdges[1] = TimeVaryingSemiring({{0,0,1,1}});
    SREdges[2] = TimeVaryingSemiring({{4,4,5,5}});
    SREdges[3] = TimeVaryingSemiring({{1,1,2,2}});
    SREdges[4] = TimeVaryingSemiring({{2,2,3,3}});

	GenericFloydWarshall<TimeVaryingSemiring> aptv(G, SREdges);


    std::cout << "HIER\n";
    aptv.run();
    std::cout << aptv.getDistance(2,1) << "\n";
    // EXPECT_EQ(aptv.getDistance(0,0), TimeVaryingSemiring({{0,0,0,0}}));
    // EXPECT_EQ(aptv.getDistance(0,1), TimeVaryingSemiring({{0,0,0,0}}));
    // EXPECT_EQ(aptv.getDistance(0,2), TimeVaryingSemiring({{0,0,0,0}}));
    // EXPECT_EQ(aptv.getDistance(0,3), TimeVaryingSemiring({{0,0,0,0}}));
};

}

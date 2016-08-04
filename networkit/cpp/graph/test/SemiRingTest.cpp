/*
 * SemiRingTest.cpp
 *
 *  Created on: 21.07.2014
 *      Author: ebergamini
 */

#include "SemiRingTest.h"
#include "../DynBFS.h"
#include "../BFS.h"
#include "../DynDijkstra.h"
#include "../SemiringDijkstra.h"
#include "../FloydWarshall.h"
#include "../BellmanFord.h"
#include "../../io/METISGraphReader.h"
#include "../../auxiliary/Log.h"

#include <stack>

#include "../ShortestPathSemiring.h"
#include "../SafestPathSemiring.h"
#include "../TimeVaryingSemiring.h"


namespace NetworKit {

TEST_F(SemiRingTest, testDijkstra) {

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
    SREdges[5] = 0.1;

	SemiringDijkstra<SafestPathSemiring> sssp(G, 0, SREdges);
    FloydWarshall<SafestPathSemiring> apapp(G, SREdges);

    sssp.run();
	apapp.run();
    std::vector<std::vector<SafestPathSemiring>> dist = apapp.getDistances();
    for (u_int i = 0; i < G.upperNodeIdBound(); i++)
    {
        for (u_int j = 0; j < G.upperNodeIdBound(); j++) {
            std::cout << " " << dist[i][j];
        }
        std::cout << "\n";
    }
    
};

}

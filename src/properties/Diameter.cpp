/*
 * Diameter.cpp
 *
 *  Created on: 21.11.13
 *      Author: lbarth, dweiss
 */

#include "Diameter.h"

namespace GrauBart {

using namespace NetworKit;

Diameter::Diameter()
{
	// TODO Auto-generated constructor stub

}

Diameter::~Diameter()
{
	// TODO Auto-generated destructor stub
}

count
Diameter::randomLowerBound(const Graph &G) 
{
	// TODO implement me
}

count
Diameter::randomUpperBound(const Graph &G) 
{
	// First, build BFS from a random node 

	// TODO make sure this is a node with high degree

	node root = random() % (G.nodeCount());

	std::vector<count> BFSdist = BFS().run(G, root);

	// Now determine vertex with maximum eccentricity in the spanning tree
	node max = 0;
	count maxDist = 0;
	node cur = 0;
	for (auto it = BFSdist.begin(); it != BFSdist.end(); ++it) {
		count dist = *it;

		if (dist > maxDist) {
			maxDist = dist;
			max = cur;
		}
		cur++;
	}

	// OK, now let's find the vertex with maximum distance vom max within the tree. 
}

} /* namespace NetworKit */

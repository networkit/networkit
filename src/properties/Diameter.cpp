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

std::pair<count,count>
estimateDiameterRange(const Graph& G) const
{
  count lower = 0;
  count upper = std::numeric_limits<count>::max();
  node root = 0;
  count dist = 0;
  std::vector<count> BFSdist;
  
  for (int i = 0; (i < k) && (upper - lower > 2); i++ {
    // improve lower bound
    root = random() % (G.nodeCount());
    BFSdist = BFS().run(G, root);
    dist = *(std::max_element(BFSdist.begin(), BFSdist.end()));
    
    if (dist > lower) {
      lower = dist;
    }
    
    // improve upper bound
    root = random() % (G.nodeCount());
    dist = getUpperBound(G, root);
    
    if (dist < upper) {
      upper = dist;
    }
  }
  
  return std::pair<count, count>(lower, upper);
}

count
getUpperBound(const Graph &G, const node root) const
{
  // TODO: build BFS tree for high-degree node
  // in BFS tree find node v with d(root,v) = ecc(root)
  // in same tree find ecc(v) and return
}

} /* namespace NetworKit */

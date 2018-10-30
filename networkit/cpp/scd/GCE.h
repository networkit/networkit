/* GCE.h
 *
 * Created on: 06.05.2013
 * Author: cls
 */


#ifndef GCE_H_
#define GCE_H_

#include <unordered_set>

#include "SelectiveCommunityDetector.h"
#include "../auxiliary/SetIntersector.h"


namespace NetworKit {


/**
 * The Greedy Community Expansion algorithm.
 *
 * Greedily adds nodes from the shell to improve community quality.
 */
class GCE: public SelectiveCommunityDetector {

public:

	GCE(const Graph& G, std::string objective);


	std::map<node, std::set<node> >  run(const std::set<node>& seeds) override;

	/**
	 * @param[in]	s	seed node
	 *
	 * @param[out]		community as a set of nodes
	 */
	std::set<node> expandSeed(node s);

protected:

    std::string objective;    // name of objective function


};

} /* namespace NetworKit */
#endif

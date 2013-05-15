/*
 * SelectiveCommunityDetector.h
 *
 *  Created on: 15.05.2013
 *      Author: cls
 */

#ifndef SELECTIVECOMMUNITYDETECTOR_H_
#define SELECTIVECOMMUNITYDETECTOR_H_

namespace NetworKit {

class SelectiveCommunityDetector {

public:

	SelectiveCommunityDetector();

	virtual ~SelectiveCommunityDetector();

	/**
	 * @param[in]	G		the graph
	 * @param[in]	seeds	set of seed nodes
	 *
	 * @param[out]			partial clustering, i.e. some entries are undefined
	 */
	virtual Clustering& run(Graph& G, std::vector<node> seeds) = 0;
};

} /* namespace NetworKit */
#endif /* SELECTIVECOMMUNITYDETECTOR_H_ */

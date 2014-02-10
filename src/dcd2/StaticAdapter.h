/*
 * StaticAdapter.h
 *
 *  Created on: 10.01.2014
 *      Author: cls
 */

#ifndef STATICADAPTER_H_
#define STATICADAPTER_H_

#include "DynCommunityDetector.h"
#include "../community/Clusterer.h"

namespace NetworKit {

/**
 * This class allows us to pass a static community detection algorithm
 * as a dynamic one.
 */
class StaticAdapter: public NetworKit::DynCommunityDetector {

public:

	StaticAdapter(Clusterer* algo);

	void update(std::vector<GraphEvent>& stream) override;

	Partition detect(bool restart) override;


private:

	Clusterer* algo;
};

} /* namespace NetworKit */

#endif /* STATICADAPTER_H_ */

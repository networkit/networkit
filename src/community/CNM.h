/*
 * CNM.h
 *
 *  Created on: Jun 10, 2013
 *      Author: michi
 */

#ifndef CNM_H_
#define CNM_H_

#include "Clusterer.h"

namespace NetworKit {

/**
 * Clustering algorithm due to Clauset, Newman and Moore.
 */
class CNM : public NetworKit::Clusterer {
public:
	CNM();
	virtual ~CNM();


	virtual Clustering run(Graph &graph);

	virtual std::string toString() {
		return "CNM";
	}
};

}

#endif /* CNM_H_ */

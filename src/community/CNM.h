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


	Clustering run(Graph &graph) override;

	std::string toString() const override {
		return "CNM";
	}
};

}

#endif /* CNM_H_ */

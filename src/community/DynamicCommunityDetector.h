/*
 * DynamicClusterer.h
 *
 *  Created on: 10.04.2013
 *      Author: cls
 */

#ifndef DYNAMICCOMMUNITYDETECTOR_H_
#define DYNAMICCOMMUNITYDETECTOR_H_

#include "../clustering/Clustering.h"
#include "../dynamics/GraphEventHandler.h"

namespace NetworKit {

class DynamicCommunityDetector: public NetworKit::GraphEventHandler {


public:

	DynamicCommunityDetector(); // nullary constructor needed for Python shell - do not create instances with this

	DynamicCommunityDetector(Graph& G);

	virtual ~DynamicCommunityDetector();

	virtual Clustering run() = 0;

	virtual std::string toString() const = 0;


protected:

	Graph* G;

};

} /* namespace NetworKit */
#endif /* DYNAMICCOMMUNITYDETECTOR_H_ */

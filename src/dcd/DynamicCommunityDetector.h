/*
 * DynamicClusterer.h
 *
 *  Created on: 10.04.2013
 *      Author: cls
 */

#ifndef DYNAMICCOMMUNITYDETECTOR_H_
#define DYNAMICCOMMUNITYDETECTOR_H_

#include <sstream>

#include "../clustering/Clustering.h"
#include "../dynamics/GraphEventHandler.h"

namespace NetworKit {

class DynamicCommunityDetector: public NetworKit::GraphEventHandler {


public:


	DynamicCommunityDetector();

	virtual ~DynamicCommunityDetector();

	/**
	 * Set the Graph instance. Needs to be called before calling run().
	 */
	virtual void setGraph(Graph& G) = 0;

	virtual Clustering run() = 0;

	virtual std::string toString() const = 0;



};

} /* namespace NetworKit */
#endif /* DYNAMICCOMMUNITYDETECTOR_H_ */

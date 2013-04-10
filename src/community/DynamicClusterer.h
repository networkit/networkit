/*
 * DynamicClusterer.h
 *
 *  Created on: 10.04.2013
 *      Author: cls
 */

#ifndef DYNAMICCLUSTERER_H_
#define DYNAMICCLUSTERER_H_

#include "../clustering/Clustering.h"
#include "../dynamics/GraphEventHandler.h"

namespace NetworKit {

class DynamicClusterer: public NetworKit::GraphEventHandler {


public:

	DynamicClusterer(Graph& G);

	virtual ~DynamicClusterer();

	virtual Clustering run();

	virtual std::string toString() const = 0;


protected:

	Graph* G;

};

} /* namespace NetworKit */
#endif /* DYNAMICCLUSTERER_H_ */

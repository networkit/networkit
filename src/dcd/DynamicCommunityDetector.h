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
#include "../auxiliary/Timer.h"

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

	virtual std::vector<count> getTimerHistory();


	/** Do not forget to also implement the event handler abstract methods from GraphEventHandler **/



protected:

	Graph* G;							//!< the current graph instance

	Aux::Timer runtime;					//!< runtime measurement
	std::vector<count> timerHistory;	//1< stores a history of runtimes



};

} /* namespace NetworKit */
#endif /* DYNAMICCOMMUNITYDETECTOR_H_ */

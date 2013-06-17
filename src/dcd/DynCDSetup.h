/*
 * DynCDSetup.h
 *
 *  Created on: 16.06.2013
 *      Author: cls
 */

#ifndef DYNCDSETUP_H_
#define DYNCDSETUP_H_

#include "DynamicCommunityDetector.h"
#include "../generators/DynamicGraphGenerator.h"

namespace NetworKit {

class DynCDSetup {

public:
	/**
	 * @param[in]	dynGen		dynamic graph generator
	 * @param[in]	dynDetectors	collection of dynamic community detection algorithms
	 * @param[in]	tMax			maximum number of time steps
	 * @param[in]	deltaT			time steps between two algorithm runs
	*/

	DynCDSetup(DynamicGraphSource& dynGen, std::vector<DynamicCommunityDetector*>& dynDetectors, count tMax, count deltaT=1);

	virtual ~DynCDSetup();

	virtual void run();

private:
	DynamicGraphSource* gen;	//!< pointer to dynamic graph generator
	std::vector<DynamicCommunityDetector*> detectors;	//!< pointer to collection of dynamic community detection algorithms
	Graph* G;
	GraphEventProxy* Gproxy;
	count deltaT;	//!< number of time steps between two algorithm runs
	count tMax;		//!< maximum number of time steps

};

} /* namespace NetworKit */
#endif /* DYNCDSETUP_H_ */

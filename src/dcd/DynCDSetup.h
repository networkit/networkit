/*
 * DynCDSetup.h
 *
 *  Created on: 16.06.2013
 *      Author: cls
 */

#ifndef DYNCDSETUP_H_
#define DYNCDSETUP_H_

#include "DynamicCommunityDetector.h"
#include "../generators/DynamicGraphSource.h"
#include "../auxiliary/Debug.h"
#include "../community/Clusterer.h"
#include "../clustering/Modularity.h"
#include "../clustering/DynamicNMIDistance.h"

namespace NetworKit {

class DynCDSetup {

public:
	/**
	 * Construct a setup.
	 *
	 * @param[in]	dynGen		dynamic graph generator
	 * @param[in]	dynDetectors	collection of dynamic community detection algorithms
	 * @param[in]	tMax			maximum number of time steps
	 * @param[in]	deltaT			time steps between two algorithm runs
	*/

	DynCDSetup(DynamicGraphSource& dynGen, std::vector<DynamicCommunityDetector*>& dynDetectors, count tMax, count deltaT=1);

	virtual ~DynCDSetup();


	virtual void setStatic(Clusterer* staticAlgo);


	/**
	 * Run the setup.
	 */
	virtual void run();

	/**
	 * Return pointer to the graph instance.
	 */
	virtual Graph* getGraph();


	/**
	 * Return copy of the graph instance.
	 */
	virtual Graph getGraphCopy();


	/**
	 * After calling this, the setup calculates modularity values.
	 */
	virtual void checkModularity();

	/**
	 * After calling this, the setup checks the number of clusters
	 */
	virtual void checkNumberOfCommunities();


	virtual void checkNMIDistance();




public:
	DynamicGraphSource* gen;	//!< pointer to dynamic graph generator
	std::vector<DynamicCommunityDetector*> detectors;	//!< pointer to collection of dynamic community detection algorithms
	Graph* G;
	GraphEventProxy* Gproxy;
	count deltaT;	//!< number of time steps between two algorithm runs
	count tMax;		//!< maximum number of time steps

	Clusterer* staticAlgo; //!< a static clustering algorithm

	std::vector<std::vector<Clustering> > results; //!< the resulting communities per algorithm per run
	std::vector<Clustering> staticClusterings; //!< if there is a static algorithm, store its results here

protected:

	bool checkMod = false; //!< if this is true, we check modularity
	bool checkNumCom = false; //!< if this is true, we check the number of communities
	bool checkNMID = false; 	//!< if this is true, we check NMIDistance between consecutive clusterings

};

} /* namespace NetworKit */
#endif /* DYNCDSETUP_H_ */

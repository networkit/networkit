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
#include "../community/Clusterer.h"


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


	virtual void checkContinuity();




public:
	DynamicGraphSource* gen;	//!< pointer to dynamic graph generator
	std::vector<DynamicCommunityDetector*> detectors;	//!< pointer to collection of dynamic community detection algorithms
	Graph* G;
	GraphEventProxy* Gproxy;
	count nDetectors; //!< number of dynamic detectors
	count deltaT;	//!< number of time steps between two algorithm runs
	count tMax;		//!< maximum number of time steps

	Clusterer* staticAlgo; //!< a static clustering algorithm

	// TIMELINES


	std::vector<count> nTimeline;
	std::vector<count> mTimeline;

	// for static
	std::vector<Clustering> staticClusteringTimeline; //!< if there is a static algorithm, store its results here
	std::vector<count> staticTimerTimeline;
	std::vector<double> staticQualityTimeline;
	std::vector<count> staticNCommunitiesTimeline;
	std::vector<double> staticContinuityTimeline;

	// once per dynamic detector
	std::vector<std::vector<Clustering> > dynamicClusteringTimelines; //!< the resulting communities per algorithm per run
	std::vector<std::vector<double> > qualityTimelines;
	std::vector<std::vector<count> > nCommunitiesTimelines;
	std::vector<std::vector<double> > continuityTimelines;

	bool checkMod; //!< if this is true, we check modularity
	bool checkNumCom; //!< if this is true, we check the number of communities
	bool checkSampledRand;
	bool checkNMID;

};

} /* namespace NetworKit */
#endif /* DYNCDSETUP_H_ */

/*
 * DynamicCommunityDetection.h
 *
 *  Created on: 10.01.2014
 *      Author: cls
 */

#ifndef DYNAMICCOMMUNITYDETECTION_H_
#define DYNAMICCOMMUNITYDETECTION_H_


#include <string>

#include "../graph/Graph.h"
 #include "../clustering/Clustering.h"

namespace NetworKit {

/*
 *
 */
class DynamicCommunityDetection {
public:
	DynamicCommunityDetection(std::string inputPath, std::string algoName, count interval=1, 
		bool recordQuality=true, bool recordContinuity=true);

	void run();

	std::vector<count> getUpdateTimeline();

	std::vector<count> getDetectTimeline();

	std::vector<double> getQualityTimeline();

	std::vector<double> getContinuityTimeline();

	std::vector<std::pair<count, count> > getGraphSizeTimeline();

private:

	std::string inputPath;
	std::string algoName;
	count interval;

	bool recordQuality;
	bool recordContinuity;

	std::vector<double> quality;
	std::vector<count> updateTime;
	std::vector<count> detectTime;
	std::vector<double> continuity;
	std::vector<std::pair<count, count> > size; // records graph size

	Clustering previous; // communities from the previous run

};

} /* namespace NetworKit */

#endif /* DYNAMICCOMMUNITYDETECTION_H_ */

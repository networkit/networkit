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
	DynamicCommunityDetection(std::string inputPath, std::string algoName, std::string updateStrategy="isolate", count interval=1,
		std::vector<std::string> recordSettings = {"quality", "graphSize"});

	void run();

	std::vector<count> getUpdateTimeline();

	std::vector<count> getDetectTimeline();

	std::vector<double> getQualityTimeline();

	std::vector<double> getContinuityTimeline();

	std::vector<std::pair<count, count> > getGraphSizeTimeline();

private:

	std::string inputPath;
	std::string algoName;
	std::string updateStrategy;
	count interval;

	std::vector<std::string> recordSettings; // determines which properties to record

	std::vector<double> quality;
	std::vector<count> updateTime;
	std::vector<count> detectTime;
	std::vector<double> continuity;
	std::vector<std::pair<count, count> > size; // records graph size
	std::vector<std::vector<count>> communitySizes;

	Clustering previous; // communities from the previous run

};

} /* namespace NetworKit */

#endif /* DYNAMICCOMMUNITYDETECTION_H_ */

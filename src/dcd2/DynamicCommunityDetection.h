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
		std::vector<std::string> recordSettings = {"quality", "continuity"});

	void run();

	std::vector<std::pair<count, count> > getGraphSizeTimeline();

	std::vector<double> getTimeline(std::string key);

private:

	std::string inputPath;
	std::string algoName;
	std::string updateStrategy;
	count interval;

	std::vector<std::string> recordSettings; // determines which properties to record

	std::vector<double> quality;
	std::vector<double> updateTime;
	std::vector<double> detectTime;
	std::vector<double> continuity;

	std::vector<std::vector<count>> communitySizes;
	std::vector<double> communityCount;

	std::vector<std::pair<count, count> > size; // records graph size


	Clustering previous; // communities from the previous run

};

} /* namespace NetworKit */

#endif /* DYNAMICCOMMUNITYDETECTION_H_ */

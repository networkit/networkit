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

namespace NetworKit {

/*
 *
 */
class DynamicCommunityDetection {
public:
	DynamicCommunityDetection(std::string inputPath, std::string algoName, count interval=1);

	void run();

	std::vector<count> getUpdateTimeline();

	std::vector<count> getDetectTimeline();

	std::vector<double> getQualityTimeline();

	std::vector<double> getContinuityTimeline();

private:

	std::string inputPath;
	std::string algoName;
	count interval;

	std::vector<double> quality;
	std::vector<count> updateTime;
	std::vector<count> detectTime;
	std::vector<double> continuity;

};

} /* namespace NetworKit */

#endif /* DYNAMICCOMMUNITYDETECTION_H_ */

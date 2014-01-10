/*
 * DynamicCommunityDetection.h
 *
 *  Created on: 10.01.2014
 *      Author: cls
 */

#ifndef DYNAMICCOMMUNITYDETECTION_H_
#define DYNAMICCOMMUNITYDETECTION_H_


#include <string>

namespace NetworKit {

/*
 *
 */
class DynamicCommunityDetection {
public:
	DynamicCommunityDetection();

	void run(std::string inputPath, std::string algoName);

};

} /* namespace NetworKit */

#endif /* DYNAMICCOMMUNITYDETECTION_H_ */

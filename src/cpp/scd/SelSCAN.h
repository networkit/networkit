/**
 * Created on: 06.05.2013
 * Author: cls
 */

#ifndef SELSCAN_H_
#define SELSCAN_H_

#include <unordered_set>

#include "SelectiveCommunityDetector.h"
#include "../distmeasures/NodeDistance.h"


namespace NetworKit {


/**
 * TODO:
 */
class SelSCAN: public NetworKit::SelectiveCommunityDetector {

public:

	SelSCAN(const Graph& G, std::string nodeDistanceMeasure, count kappa, double epsilon);


	void run(std::set<unsigned int>& seeds) override;

protected:

	NodeDistance* nodeDist;
	count kappa;
	double epsilon;


};

} /* namespace NetworKit */
#endif
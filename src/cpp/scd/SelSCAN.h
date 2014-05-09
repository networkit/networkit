/**
 * Created on: 06.05.2013
 * Author: cls
 */

#ifndef SELSCAN_H_
#define SELSCAN_H_

#include <unordered_set>

#include "SelectiveCommunityDetector.h"
#include "../distmeasures/AlgebraicDistance.h"


namespace NetworKit {


/**
 * TODO:
 */
class SelSCAN: public NetworKit::SelectiveCommunityDetector {

public:

	SelSCAN(const Graph& G, count kappa, double epsilon);

	SelSCAN(const Graph& G, count kappa, double epsilon, AlgebraicDistance& ad);


	void run(std::set<unsigned int>& seeds) override;

protected:

	count kappa;
	double epsilon;
	AlgebraicDistance* algebraicDistance;
	std::function<double(node, node)> dist; 	// node distance measure


};

} /* namespace NetworKit */
#endif

/*
 * DynPLM.h
 *
 *  Created on: 03.01.2014
 *      Author: cls
 */

#ifndef DYNPLM_H_
#define DYNPLM_H_

#include "DynCommunityDetector.h"

#include "../community/PLM2.h"

namespace NetworKit {

class DynPLM: public NetworKit::DynCommunityDetector {

public:

	/**
	 * @param[in]	refine	add a second move phase to refine the communities
	 * @param[in]	par		parallelization strategy
	 * @param[in]	gamma	multi-resolution modularity parameter:
	 * 							1.0 -> standard modularity
	 * 							0.0 -> one community
	 * 							2m 	-> singleton communities
	 *
	 */
	DynPLM(std::string prepStrategy="isolate", bool refine=false, double gamma = 1.0, std::string par="balanced");

	void update(std::vector<GraphEvent>& stream) override;

	Clustering detect() override;

private:

	Clustering run(Graph& G);

	std::string prepStrategy; // determines which dynamization strategy to use
	std::string parallelism;
	bool refine;
	double gamma = 1.0;

	PLM2 plm2; // instance of the static algorithm

};

} /* namespace NetworKit */

#endif /* DYNPLM_H_ */

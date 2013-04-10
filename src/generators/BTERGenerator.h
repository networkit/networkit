/*
 * BTERGenerator.h
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#ifndef BTERGENERATOR_H_
#define BTERGENERATOR_H_

#include "StaticGraphGenerator.h"

namespace NetworKit {

/**
 * Block Two-level Erdos-Renyi Graph Model
 *
 * Paper: http://www.sandia.gov/~tgkolda/pubs/bibtgkfiles/BTER-arXiv-1302.6636v1.pdf
 * Video introduction: http://www.youtube.com/watch?v=kF-pKfg846M
 */
class BTERGenerator: public NetworKit::StaticGraphGenerator {

public:

	/**
	 * @param[in]	degreeDistribution			number of nodes per degree
	 * @param[in]	clusteringCoefficients		mean clusering coefficient for nodes per degree
	 */
	BTERGenerator(std::vector<count> degreeDistribution, std::vector<count> clusteringCoefficients);

	virtual ~BTERGenerator();

	virtual Graph generate();


protected:

	virtual void preprocessing();

	virtual void phaseOne();

	virtual void phaseTwo();

};


} /* namespace NetworKit */
#endif /* BTERGENERATOR_H_ */

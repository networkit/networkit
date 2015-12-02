/*
 * ChungLuScore.h
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#ifndef CHUNGLUSCORE_H
#define CHUNGLUSCORE_H

#include "../edgescores/EdgeScore.h"

namespace NetworKit {

class ChungLuScore : public EdgeScore<double> {

public:
	ChungLuScore(const Graph& graph);
	
	virtual void run() override;
};

} /* namespace NetworKit */

#endif // CHUNGLUSCORE_H

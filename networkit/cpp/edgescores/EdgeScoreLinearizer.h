/*
 * EdgeScoreLinearizer.h
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#ifndef EDGESCORELINEARIZER_H
#define EDGESCORELINEARIZER_H

#include "EdgeScore.h"

namespace NetworKit {

class EdgeScoreLinearizer : public EdgeScore<double> {

private:
	const std::vector<double>& attribute;
	bool inverse;

public:
	EdgeScoreLinearizer(const Graph& graph, const std::vector<double>& attribute, bool inverse = false);

	virtual double score(edgeid eid) override;
	virtual double score(node u, node v) override;
	virtual void run() override;

};

} // namespace NetworKit

#endif // EDGESCORELINEARIZER_H

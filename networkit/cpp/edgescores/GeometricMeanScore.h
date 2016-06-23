/*
 * GeometricMeanScore.h
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#ifndef GEOMETRICMEANSCORE_H
#define GEOMETRICMEANSCORE_H

#include "EdgeScore.h"

namespace NetworKit {

class GeometricMeanScore : public EdgeScore<double> {

private:
	const std::vector<double>& attribute;

public:
	GeometricMeanScore(const Graph& G, const std::vector<double>& attribute);
	virtual double score(edgeid eid) override;
	virtual double score(node u, node v) override;
	virtual void run() override;
};

} // namespace NetworKit

#endif // GEOMETRICMEANSCORE_H

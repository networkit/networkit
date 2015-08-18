/*
 * LocalDegreeScore.h
 *
 *  Created on: 28.08.2014
 *      Author: Gerd Lindner
 */

#ifndef LOCALDEGREESCORE_H_
#define LOCALDEGREESCORE_H_

#include "../edgescores/EdgeScore.h"

namespace NetworKit {

/**
 * Local Degree sparsification method.
 */
class LocalDegreeScore : public EdgeScore<double> {

public:

	LocalDegreeScore(const Graph& graph);
	virtual std::vector<double> scores() override;
	virtual double score(edgeid eid) override;
	virtual double score(node u, node v) override;
	virtual void run() override;

private:
	const Graph& graph;
	std::vector<double> scoreData;

};

}
/* namespace NetworKit */
#endif /* LOCALDEGREEATTRIBUTIZER_H_ */

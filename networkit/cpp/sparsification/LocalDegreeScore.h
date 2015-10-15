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
 * See 'Structure-Preserving Sparsification of Social Networks' by Lindner, Staudt, Hamann.
 */
class LocalDegreeScore : public EdgeScore<double> {

public:

	LocalDegreeScore(const Graph& G);
	virtual void run() override;
	virtual double score(edgeid eid) override;
	virtual double score(node u, node v) override;

};

}
/* namespace NetworKit */
#endif /* LOCALDEGREESCORE_H_ */

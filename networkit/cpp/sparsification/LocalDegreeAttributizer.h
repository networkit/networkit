/*
 * LocalDegreeAttributizer.h
 *
 *  Created on: 28.08.2014
 *      Author: Gerd Lindner
 */

#ifndef LOCALDEGREEATTRIBUTIZER_H_
#define LOCALDEGREEATTRIBUTIZER_H_

#include "../edgescores/EdgeAttribute.h"

namespace NetworKit {

/**
 * EXPERIMENTAL
 */
class LocalDegreeAttributizer : public EdgeAttribute<double> {

public:

	LocalDegreeAttributizer(const Graph& graph);
	virtual std::vector<double> getAttribute() override;

private:
	const Graph& graph;

};

}
/* namespace NetworKit */
#endif /* LOCALDEGREEATTRIBUTIZER_H_ */

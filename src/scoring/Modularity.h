/*
 * Modularity.h
 *
 *  Created on: 15.10.2012
 *      Author: cls
 */

#ifndef MODULARITY_H_
#define MODULARITY_H_

#include "EdgeScoring.h"

namespace EnsembleClustering {

class Modularity: public EnsembleClustering::EdgeScoring {

public:

	Modularity();

	virtual ~Modularity();

	virtual double scoreEdge(id u, id v);
};

} /* namespace EnsembleClustering */
#endif /* MODULARITY_H_ */

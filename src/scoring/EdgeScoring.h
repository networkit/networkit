/*
 * EdgeScoring.h
 *
 *  Created on: 15.10.2012
 *      Author: cls
 */

#ifndef EDGESCORING_H_
#define EDGESCORING_H_

namespace EnsembleClustering {

class EdgeScoring {

public:

	EdgeScoring();

	virtual ~EdgeScoring();

	virtual double scoreEdge(id u, id v);
};

} /* namespace EnsembleClustering */
#endif /* EDGESCORING_H_ */

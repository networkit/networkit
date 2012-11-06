/*
 * StaticGraphImplementor.h
 *
 *  Created on: 06.11.2012
 *      Author: cls
 */

#ifndef STATICGRAPHIMPLEMENTOR_H_
#define STATICGRAPHIMPLEMENTOR_H_

#include "GraphImplementor.h"

namespace EnsembleClustering {

class StaticGraphImplementor: public EnsembleClustering::GraphImplementor {
public:
	StaticGraphImplementor();
	virtual ~StaticGraphImplementor();
};

} /* namespace EnsembleClustering */
#endif /* STATICGRAPHIMPLEMENTOR_H_ */

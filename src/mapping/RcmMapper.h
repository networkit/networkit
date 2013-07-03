/*
 * RcmMapper.h
 *
 *  Created on: Jul 2, 2013
 *      Author: Henning
 */

#ifndef RCMMAPPER_H_
#define RCMMAPPER_H_

#include "StaticMapper.h"

namespace NetworKit {

class RcmMapper: public NetworKit::StaticMapper {
public:
	RcmMapper();
	virtual ~RcmMapper();
};

} /* namespace NetworKit */
#endif /* RCMMAPPER_H_ */

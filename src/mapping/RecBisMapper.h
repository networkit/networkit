/*
 * RecBisMapper.h
 *
 *  Created on: Jul 2, 2013
 *      Author: Henning
 */

#ifndef RECBISMAPPER_H_
#define RECBISMAPPER_H_

#include "StaticMapper.h"

namespace NetworKit {

class RecBisMapper: public NetworKit::StaticMapper {
public:
	RecBisMapper();
	virtual ~RecBisMapper();
};

} /* namespace NetworKit */
#endif /* RECBISMAPPER_H_ */

/*
 * METIStoSTINGER.h
 *
 *  Created on: 27.11.2012
 *      Author: cls
 */

#ifndef METISTOSTINGER_H_
#define METISTOSTINGER_H_

#include <string>

extern "C" {
	#include "stinger.h"
}

typedef stinger graph;


namespace EnsembleClustering {

class METIStoSTINGER {

public:

	METIStoSTINGER();

	virtual ~METIStoSTINGER();

	virtual graph* read(std::string graphPath);
};

} /* namespace EnsembleClustering */
#endif /* METISTOSTINGER_H_ */

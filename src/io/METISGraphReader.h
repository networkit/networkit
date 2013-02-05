/*
 * METISGraphReader.h
 *
 *  Created on: 17.01.2013
 *      Author: cls
 */

#ifndef METISGRAPHREADER_H_
#define METISGRAPHREADER_H_

#include "GraphReader.h"

#include "METISParser.h"

namespace EnsembleClustering {

class METISGraphReader: public EnsembleClustering::GraphReader {

public:

	METISGraphReader();

	virtual ~METISGraphReader();

	virtual Graph read(std::string path);
};

} /* namespace EnsembleClustering */
#endif /* METISGRAPHREADER_H_ */

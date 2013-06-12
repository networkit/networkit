/*
 * DibapGraphReader.h
 *
 *  Created on: Jun 12, 2013
 *      Author: Henning
 */

#ifndef DIBAPGRAPHREADER_H_
#define DIBAPGRAPHREADER_H_

#include "GraphReader.h"
#include "../graph/Graph.h"
#include <cstdio>


// codes in file headers to distinguish type
#define IO_TYPE_XX (('X' << 8) | 'X')
#define IO_TYPE_GI (('G' << 8) | 'I')
#define IO_TYPE_GF (('G' << 8) | 'F')
#define IO_TYPE_HI (('H' << 8) | 'I')
#define IO_TYPE_HF (('H' << 8) | 'F')
#define IO_TYPE_P2 (('P' << 8) | '2')
#define IO_TYPE_P4 (('P' << 8) | '4')
#define IO_TYPE_AA (('A' << 8) | 'A')
#define IO_TYPE_T2 (('T' << 8) | '2')
#define IO_TYPE_TE (('T' << 8) | 'E')


namespace NetworKit {

class DibapGraphReader: public NetworKit::GraphReader {
public:
	DibapGraphReader();
	virtual ~DibapGraphReader();

	virtual Graph read(std::string path);
};

} /* namespace NetworKit */
#endif /* DIBAPGRAPHREADER_H_ */

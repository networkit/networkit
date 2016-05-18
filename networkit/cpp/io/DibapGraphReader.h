/*
 * DibapGraphReader.h
 *
 *  Created on: Jun 12, 2013
 *      Author: Henning
 */

#ifndef DIBAPGRAPHREADER_H_
#define DIBAPGRAPHREADER_H_

#if !defined _WIN32 && !defined _WIN64 && !defined WIN32 && !defined WIN64

#include "GraphReader.h"
#include "../graph/Graph.h"

#include <cstdio>
#include <netinet/in.h>

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

/**
 * @ingroup io
 * TODO: class documentation
 */
class DibapGraphReader: public NetworKit::GraphReader {
public:
	DibapGraphReader() = default;

	virtual Graph read(const std::string& path) override;
};

} /* namespace NetworKit */

#endif /* check for non-Windows */

#endif /* DIBAPGRAPHREADER_H_ */

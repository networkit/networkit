/*
 * METISGraphWriter.h
 *
 *  Created on: 30.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef GMLGRAPHWRITER_H_
#define GMLGRAPHWRITER_H_

#include <fstream>

#include "GraphWriter.h"

namespace NetworKit {

/**
 * @ingroup io
 * Writes a (so far unweighted) graph and its coordinates as a GML file.
 * TODO: write unit test
 */
class GMLGraphWriter: public NetworKit::GraphWriter {
public:
	/** Default constructor */
	GMLGraphWriter() = default;

	/**
	 * Write a graph @a G and its coordinates to a GML file.
	 *
	 * @param[in]	G		Graph of type NetworKit with 2D coordinates.
	 * @param[in]	path	Path to file.
	 */
	virtual void write(Graph& G, const std::string& path) override;
};

} /* namespace NetworKit */
#endif /* GMLGRAPHWRITER_H_ */

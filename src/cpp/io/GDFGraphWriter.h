/*
 * GDFGraphWriter.h
 *
 *  Created on: 27.10.2013
 *      Author: Stefan Bertsch
 */

#ifndef GDFGRAPHWRITER_H_
#define GDFGRAPHWRITER_H_

#include <fstream>

#include "GraphWriter.h"

namespace NetworKit {

/**
 * @ingroup io
 * Writes a graph and its coordinates to a file in GDF format.[1]
 * [1] http://guess.wikispot.org/The_GUESS_.gdf_format
 * TODO: write unit test
 */
class GDFGraphWriter: public NetworKit::GraphWriter {
protected:
	/**
	 * Write a graph @a G and its coordinates to a file of GDF format.
	 *
	 * @param[in]	G			Graph of type Networkit with 2D or 3D coordinates.
	 * @param[in]	weighted	@c true if @a G is weighted (makes no difference whether is true or not)
	 * @param[in]	path		Path to file.
	 * @param[in]	dim			Dimension of coordinates.
	 */
	virtual void writeGeneric(Graph& G, bool weighted, const std::string& path, count dim);


public:
	/** Default destructor */
	GDFGraphWriter() = default;
	
	/**
	 * Write a graph @a G and its coordinates to a file of GDF format.
	 *
	 * @param[in]	G		Graph of type NetworKit with 2D coordinates.
	 * @param[in]	path	Path to file.
	 */
	virtual void write(Graph& G, const std::string& path) override;
	
	/**
	 * Write a graph @a G and its coordinates to a file of GDF format.
	 *
	 * @param[in]	G			Graph of type NetworKit with 2D coordinates.
	 * @param[in]	weighted	@c true if @a G is weighted (makes no difference whether is true or not).
	 * @param[in]	path		Path to file.
	 */
	virtual void write(Graph& G, bool weighted, const std::string& path);

	/**
	 * Write a graph @a G and its coordinates to a file of GDF format.
	 *
	 * @param[in]	G		Graph of type NetworKit with 3D coordinates.
	 * @param[in]	weighted	@c true if @a G is weighted (makes no difference whether is true or not).
	 * @param[in]	path		Path to file.
	 */
	virtual void write3D(Graph& G, bool weighted, const std::string& path);
};

} /* namespace NetworKit */
#endif /* GDFGRAPHWRITER_H_ */

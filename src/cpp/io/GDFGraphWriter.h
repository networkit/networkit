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
 * Writes a graph and its coordinates to a file in GDF format.
 * TODO: write unit test
 */
class GDFGraphWriter: public NetworKit::GraphWriter {
protected:
	/**
	 * write a graph G and its coordinates to a file of GDF format.
	 *
	 * @param[in]	G			Graph of type Networkit with 2D or 3D coordinates
	 * @param[in]	weighted	true if G is weighted (makes no difference whether is true or not)
	 * @param[in]	path		path to file
	 * @param[in]	dim			dimension of coordinates
	 */
	virtual void writeGeneric(Graph& G, bool weighted, const std::string& path, count dim);


public:
	GDFGraphWriter() = default;
	
	/**
	 * write a graph G and its coordinates to a file of GDF format.
	 *
	 * @param[in]	G		Graph of type Networkit with 2D coordinates
	 * @param[in]	path	path to file
	 */
	virtual void write(Graph& G, const std::string& path) override;
	
	/**
	 * write a graph G and its coordinates to a file of GDF format.
	 *
	 * @param[in]	G			Graph of type Networkit with 2D coordinates
	 * @param[in]	weighted	true if the G is weighted (makes no difference wether is true or not)
	 * @param[in]	path		path to file
	 */
	virtual void write(Graph& G, bool weighted, const std::string& path);

	/**
	 * write a graph G and its coordinates to a file of GDF format.
	 *
	 * @param[in]	G		Graph of type Networkit with 3D coordinates
	 * @param[in]	weighted	true if the G is weighted (makes no difference wether is true or not)
	 * @param[in]	path		path to file
	 */
	virtual void write3D(Graph& G, bool weighted, const std::string& path);
};

} /* namespace NetworKit */
#endif /* GDFGRAPHWRITER_H_ */

/*
 * VNAGraphWriter.h
 *
 *  Created on: 23.10.2013
 *      Author: Stefan Bertsch
 */

#ifndef VNAGRAPHWRITER_H_
#define VNAGRAPHWRITER_H_

#include <fstream>

#include "GraphWriter.h"
#include "../structures/Partition.h"

namespace NetworKit {

/**
 * @ingroup io
 * Writes a graph and its coordinates as VNA file.
 * The VNA format is commonly used by Netdraw, and is very similar to Pajek format. It defines nodes and edges (ties),
 * and supports attributes. Each section of the file is separated by an asterisk. 
 */
class VNAGraphWriter: public NetworKit::GraphWriter {
protected:
	/**
	 * Write a graph @a G and its coordinates to a VNA file.
	 *
	 * @param[in]	G			Graph of type NetworKit with 2D coordinates.
	 * @param[in]	weighted	@c true if @a G is weighted (makes no difference whether is true or not).
	 * @param[in]	path		Path to file.
	 * @param[in]	partition	Partition of @a G (only used if @a dim == 0).
	 * @param[in]	dim			Dimension of coordinates.
	 */
	virtual void writeGeneric(Graph& G, bool weighted, const std::string& path,
	                          Partition& partition, count dim);


public:
	/**
	 * Write a graph @a G and its coordinates to a VNA file.
	 *
	 * @param[in]	G			Graph of type NetworKit with 3D coordinates.
	 * @param[in]	path		Path to file.
	 */
	virtual void write(Graph& G, const std::string& path) override;

	/**
	 * Write a graph @a G and its coordinates to a VNA file.
	 *
	 * @param[in]	G			Graph of type NetworKit with 2D coordinates.
	 * @param[in]	weighted	@c true if @a G is weighted (makes no difference whether is true or not).
	 * @param[in]	path		Path to file.
	 */
	virtual void write(Graph& G, bool weighted, const std::string& path);

	/**
	 * Write a graph @a G and its coordinates including node color to a VNA file.
	 *
	 * @param[in]	G			Graph of type NetworKit with 3D coordinates.
	 * @param[in]	weighted	@c true if the G is weighted (makes no difference whether is true or not).
	 * @param[in]	path		Path to file.
	 * @param[in]	partition	Proper Partition of @a G.
	 */
	virtual void write(Graph& G, bool weighted, const std::string& path,
	                   Partition& partition);

	/**
	 * Write a graph @a G and its coordinates to a VNA file.
	 *
	 * @param[in]	G			Graph of type NetworKit with 3D coordinates.
	 * @param[in]	weighted	@c true if @a G is weighted (makes no difference whether is true or not).
	 * @param[in]	path		Path to file.
	 */
	virtual void write3D(Graph& G, bool weighted, const std::string& path);
};

} /* namespace NetworKit */
#endif /* VNAGRAPHWRITER_H_ */

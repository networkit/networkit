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
 * Writes a graph and its coordinates as VNA file.
 * The VNA format is commonly used by Netdraw, and is very similar to Pajek format. It defines nodes and edges (ties),
 * and supports attributes. Each section of the file is separated by an asterisk. 
 */
class VNAGraphWriter: public NetworKit::GraphWriter {
protected:
	/**
	 * write a graph G and its coordinates to a VNA file.
	 *
	 * @param[in]	G			Graph of type NetworKit with 2D coordinates
	 * @param[in]	weighted	true if the G is weighted (makes no difference whether is true or not)
	 * @param[in]	path		path to file
	 * @param[in]	partition	partition of G (only used if dim == 0)
	 * @param[in]	dim			dimension of coordinates
	 */
	virtual void writeGeneric(Graph& G, bool weighted, std::string path, Partition& partition, count dim);


public:

	VNAGraphWriter();

	virtual ~VNAGraphWriter();

	/**
	 * write a graph G and its coordinates to a VNA file.
	 *
	 * @param[in]	G			Graph of type Networkit with 3D coordinates
	 * @param[in]	path		path to file
	 */
	virtual void write(Graph& G, std::string path);

	/**
	 * write a graph G and its coordinates to a VNA file.
	 *
	 * @param[in]	G			Graph of type Networkit with 2D coordinates
	 * @param[in]	weighted	true if the G is weighted (makes no difference whether is true or not)
	 * @param[in]	path		path to file
	 */
	virtual void write(Graph& G, bool weighted, std::string path);

	/**
	 * write a graph G and its coordinates including node color to a VNA file.
	 *
	 * @param[in]	G			Graph of type Networkit with 3D coordinates
	 * @param[in]	weighted	true if the G is weighted (makes no difference whether is true or not)
	 * @param[in]	path		path to file
	 * @param[in]	partition	proper Partition of G
	 */
	virtual void write(Graph& G, bool weighted, std::string path, Partition& partition);

	/**
	 * write a graph G and its coordinates to a VNA file.
	 *
	 * @param[in]	G			Graph of type Networkit with 3D coordinates
	 * @param[in]	weighted	true if the G is weighted (makes no difference whether is true or not)
	 * @param[in]	path		path to file
	 */
	virtual void write3D(Graph& G, bool weighted, std::string path);
};

} /* namespace NetworKit */
#endif /* VNAGRAPHWRITER_H_ */

/*
 * METISGraphParser.h
 *
 *  Created on: 16.10.2012
 *      Author: cls
 */

#ifndef METISGRAPHPARSER_H_
#define METISGRAPHPARSER_H_

#include <iostream>
#include <vector>
#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"
#include "log4cxx/helpers/exception.h"

#include "Graph.h"
#include "EdgeTripleGraphData.h"

namespace EnsembleClustering {


typedef unsigned int id;


/**
 * A parser for the METIS graph file format.
 *
		**METIS graph file format**


		- ASCII
		- a graph of N nodes is stored in a file of N+1 lines;
		- the first line lists the number of nodes and the number of edges;
		- If the first line contains more than two values, the extra values indicate the weights;
		- each subsequent line lists the "neighbors" of a node;
		- comment lines begin with a "%" sign;
 *
 */
class METISGraphParser {

public:

	METISGraphParser();

	virtual ~METISGraphParser();


	/**
	 * Parses and returns a graph from a METIS graph file.
	 *
	 * @param[in]	path	path of the METIS graph input file
	 *
	 * @param[out]	graph	the graph contained in the input file
	 */
	virtual Graph parse(std::string path);

	// members
	log4cxx::LoggerPtr logger; //!< logger



protected:


	/**
	 * Initializes the graph data structure before parsing.
	 */
	virtual void initGraph(int n, int m);

	/**
	 * Connects node with to all its neighbors.
	 *
	 * @param[in]	v			node id
	 * @param[in]	indices		neighbor node ids
	 */
	virtual void connectNode(id v, std::vector<id> indices);


private:

	Graph *graph; //!< output of parsing, lightweight graph interface

	EdgeTripleGraphData *graphData; //!< the actual data structure storing the graph




};

} /* namespace EnsembleClustering */
#endif /* METISGRAPHPARSER_H_ */

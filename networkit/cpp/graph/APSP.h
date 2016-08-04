/*
 * APSP.h
 *
 *  Created on: 07.07.2015
 *      Author: Arie Slobbe
 */

#ifndef APSP_H_
#define APSP_H_

#include "APAPP.h"
#include "../auxiliary/Log.h"
#include "Dijkstra.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Class for all-pair shortest path algorithm.
 */
class APSP: public APAPP<edgeweight> {

public:

	/**
	 * Creates the APSP class for @a G.
	 *
	 * @param G The graph.
	 */
	APSP(const Graph& G) : APAPP(G) {
    }

	virtual ~APSP() = default;

	/** Computes the shortest paths from each node to all other nodes. */
	void run() override {
        std::vector<edgeweight> distanceVector(G.upperNodeIdBound(), 0.0);
        distances.resize(G.upperNodeIdBound(), distanceVector);
        G.parallelForNodes([&](node u){
            Dijkstra dijk(G, u);
            dijk.run();
            distances[u] = dijk.getDistances();
        });
        hasRun = true;
    }

	/**
	* @return string representation of algorithm and parameters.
	*/
	virtual std::string toString() const override {
        return "All-Pairs Shortest Path Algorithm";
    }

protected:

	using APAPP<edgeweight>::G;
	using APAPP<edgeweight>::distances;
};

} /* namespace NetworKit */

#endif /* APSP_H_ */

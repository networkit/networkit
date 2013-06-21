/*
 * CNM.h
 *
 *  Created on: Jun 10, 2013
 *      Author: michi
 */

#ifndef CNM_H_
#define CNM_H_

#include "Clusterer.h"

namespace NetworKit {

/**
 * Clustering algorithm due to Clauset, Newman and Moore.
 * TODO: Implementation requires faster data structure for handling next merge
 * (priority queue with internal access)
 */
class CNM : public NetworKit::Clusterer {
public:
	CNM();
	virtual ~CNM();


	virtual Clustering run(Graph &graph);

	virtual std::string toString() {
		return "CNM";
	}


	class Edge {
		public:
			Edge() {}
			Edge(node u, node v, double delta) : u(u), v(v), delta(delta) {}
			virtual ~Edge() {}

			node getU();
			node getV();
			double getDelta();

			bool operator>(const Edge &e1) const {
				return delta > e1.delta;
			}

		private:
			node u;
			node v;
			double delta;
	};
};

}

#endif /* CNM_H_ */

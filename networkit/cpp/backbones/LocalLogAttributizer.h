#ifndef LOCALLOGATTRIBUTIZER_H
#define LOCALLOGATTRIBUTIZER_H

#include "AttributeGenerator.h"

namespace NetworKit {

template<typename InType>
class LocalLogAttributizer : public AttributeGenerator<InType, double> {
public:
	virtual std::vector< double > getAttribute(const Graph &g, const std::vector< InType > &attribute) override {
		/*
		* For each edge, we calculate the minimum required sparsification exponent e
		* such that the edge is contained in the backbone.
		*/

		std::vector<double> sparsificationExp(g.upperEdgeIdBound(), 1.0);

		g.parallelForNodes([&](node i) {
			count d = g.degree(i);

			/*
			 * The top d^e edges (sorted by similarity in descending order)
			 * are to be kept in the backbone.
			 */

			std::vector<std::pair<double, index> > neighbors;
			g.forNeighborsOf(i, [&](node _i, node j, edgeid eid) {
				neighbors.push_back(std::make_pair(attribute[eid], eid));
			});
			std::sort(neighbors.begin(), neighbors.end(), std::greater<std::pair<double, index> >());

			count rank = 0;
			#pragma omp critical // each value is set twice, the value can be wrong if the wrong thread wins

			for (auto it : neighbors) {
				edgeid eid = it.second;

				double e = 0.0;

				if (d > 1)
					e = rank * 1.0 / (d - 1); // Keep top 1 + e * (d-1) edges

				if (e < sparsificationExp[eid]) {
					sparsificationExp[eid] = e; // do not always write in order to avoid cache synchronization
				}

				rank++;
			}

		});

		return sparsificationExp;

	}
};
} // namespace NetworKit

#endif // LOCALLOGATTRIBUTIZER_H
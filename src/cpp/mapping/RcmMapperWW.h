/*
 * RcmMapper.h
 *
 *  Created on: Jul 2, 2013
 *      Author: Henning
 */

#ifndef RCMMAPPERWW_H_
#define RCMMAPPERWW_H_

#include "StaticMapper.h"
#include <vector>

namespace NetworKit {

typedef std::vector<node> permutation;

/**
 * TODO: class documentation
 */
class RcmMapperWW: public NetworKit::StaticMapper {
public:
	RcmMapperWW();
	virtual ~RcmMapperWW();

	std::map<index, index> run(Graph& guest, Graph& host);


private:
	permutation map(Graph& g);
	permutation invert(permutation pi);
	permutation compose(permutation pi1, permutation pi2);

	class Comparator {
		public:
			Comparator(Graph& g) : g(g) {}

			bool operator()(const node& u, const node& v) {
				return g.weightedDegree(u) < g.weightedDegree(v);
			}
		private:
			Graph& g;
	};
};

} /* namespace NetworKit */
#endif /* RCMMAPPER_H_ */

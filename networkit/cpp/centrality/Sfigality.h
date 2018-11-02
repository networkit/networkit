/*
 * Sfigality.h
 *
 *  Created on: 20.01.2016
 *      Author:Christian Staudt
 */

#ifndef SFIGALITY_H_
#define SFIGALITY_H_

#include "Centrality.h"

namespace NetworKit {

/**
 * @ingroup centrality
 * A
 */
class Sfigality: public Centrality {
public:
	/**
	 * Constructs the Sfigality class for the given Graph @a G. Sfigality is a new type of
	 * node centrality measures that is high if neighboring nodes have a higher degree, e.g. in social networks, if your friends have more friends than you. Formally:
	 *
	 * $$\sigma(u) = \frac{| \{ v: \{u,v\} \in E, deg(u) < deg(v) \} |}{ deg(u) }$$
	 *
	 * @param G The graph.

	 */
	Sfigality(const Graph& G);

	void run() override;

	/**
	 * @return the theoretical maximum degree centrality, which is $n - 1$ (including the possibility of a self-loop)
	 */
	double maximum() override;
};

} /* namespace NetworKit */

#endif /* SFIGALITY_H_ */

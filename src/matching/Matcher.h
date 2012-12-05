/*
 * Matcher.h
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#ifndef MATCHER_H_
#define MATCHER_H_

namespace EnsembleClustering {

// TODO: imports
class Matching;
class Graph;


class Matcher {

public:

	Matcher();

	virtual ~Matcher();

	virtual Matching& run(const Graph& G) = 0;
};

} /* namespace EnsembleClustering */
#endif /* MATCHER_H_ */

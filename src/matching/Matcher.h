/*
 * Matcher.h
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#ifndef MATCHER_H_
#define MATCHER_H_

namespace EnsembleClustering {

class Matcher {

public:

	Matcher();

	virtual ~Matcher();

	Matching run(Graph G) = 0;
};

} /* namespace EnsembleClustering */
#endif /* MATCHER_H_ */

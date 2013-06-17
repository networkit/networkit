/*
 * GeneratorsBenchmark.h
 *
 *  Created on: May 29, 2013
 *      Author: forigem
 */

#ifndef NOGTEST

#ifndef GENERATORSBENCHMARK_H_
#define GENERATORSBENCHMARK_H_

#include <gtest/gtest.h>


#include "../StaticBarabasiAlbertGenerator.h"

namespace NetworKit {

class GeneratorsBenchmark: public testing::Test {
public:
	GeneratorsBenchmark();
	virtual ~GeneratorsBenchmark();
};

} /* namespace NetworKit */
#endif /* GENERATORSBENCHMARK_H_ */

#endif /*NOGTEST */

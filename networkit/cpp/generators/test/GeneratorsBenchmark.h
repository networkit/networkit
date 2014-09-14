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

#include "../../auxiliary/Timer.h"

namespace NetworKit {

class GeneratorsBenchmark: public testing::Test {
protected:
	template <typename L>
	uint64_t timeOnce(L f) {
		// TODO should be moved somewhere else (Benchmark parent class or the Timer class itself)
		Aux::Timer timer;
		timer.start();
		f();
		timer.stop();
		return timer.elapsedMilliseconds();
	}
};

} /* namespace NetworKit */
#endif /* GENERATORSBENCHMARK_H_ */

#endif /*NOGTEST */

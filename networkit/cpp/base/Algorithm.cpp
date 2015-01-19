#include "Algorithm.h"
#include "../auxiliary/SignalHandling.h"
#include "../auxiliary/Log.h"
#include <exception>

namespace NetworKit {
	Algorithm::Algorithm() : hasRun(false) {

	}

	void Algorithm::run() {
		Aux::SignalHandling::init((uint64_t)this);
		try {
			runImpl();
			hasRun = true;
		} catch (Aux::InterruptException e) {
			ERROR(e.what());
			hasRun = false;
		} catch (std::runtime_error e) {
			ERROR("received runtime error:\t",e.what());
			hasRun = false;
		}
		Aux::SignalHandling::reset((uint64_t)this);
	}



} /* NetworKit */
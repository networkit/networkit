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
		} catch (std::runtime_error e) {
			ERROR(e.what());
			hasRun = false;
		}
		Aux::SignalHandling::reset((uint64_t)this);
	}

	bool Algorithm::isRunning() {
		if (!Aux::SignalHandling::gotSIGINT()) {
			return true;
		} else {
			throw std::runtime_error("Received SIGINT, aborting calculations.");
		}
	}



} /* NetworKit */
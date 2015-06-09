#include "Algorithm.h"
#include "../auxiliary/SignalHandling.h"
#include "../auxiliary/Log.h"
#include <exception>

namespace NetworKit {
	Algorithm::Algorithm() : hasRun(false) {

	}

	void Algorithm::run() {
		//Aux::SignalHandling::init((uint64_t)this);
		try {
			runImpl();
			hasRun = true;
		} catch (Aux::SignalHandling::InterruptException e) {
			ERROR(toString()," at ",this,":\t",e.what());
			hasRun = false;
		} catch (std::runtime_error e) {
			ERROR(toString()," at ",this,":\t",e.what());
			hasRun = false;
		}
		//Aux::SignalHandling::reset((uint64_t)this);
	}

	void Algorithm::assureRunning() {
		if (Aux::SignalHandling::gotSIGINT()) {
			throw Aux::SignalHandling::InterruptException();
		}
	}

	bool Algorithm::hasFinished() const {
		return hasRun;
	}

	std::string Algorithm::toString() const {
		return "Algorithm base class";
	}

} /* NetworKit */
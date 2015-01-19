#include "SignalHandling.h"
#include <exception>

namespace Aux {
bool SignalHandling::running = true;
bool SignalHandling::handlerInitialized = false;
uint64_t SignalHandling::root = 0;
bool SignalHandling::rootSet = false;

bool SignalHandling::isRunning() {
	InterruptException ie;
	if (!running) throw ie;
	return running;
}

void SignalHandling::setRunning(bool isRunning) {
	running = isRunning;
}

void SignalHandling::init(uint64_t caller) {
	if (!rootSet) {
		root = caller;
		if (!handlerInitialized) {
			auto signalHandler = [](int signum ) {
				setRunning(false);
			};
			signal(SIGINT,signalHandler);
			handlerInitialized = true;
		}
		rootSet = true;
	}
}

void SignalHandling::reset(uint64_t caller) {
	if (root == caller) {
		rootSet = false;
		running = true;
		root = 0;
	}
}

} /* Aux */
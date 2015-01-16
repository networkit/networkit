#include "SignalHandling.h"

namespace Aux {
bool SignalHandling::running = true;
bool SignalHandling::handlerInitialized = false;

bool SignalHandling::isRunning() {
	return running;
}

void SignalHandling::setRunning(bool isRunning) {
	running = isRunning;
}

void SignalHandling::init() {
	if (!handlerInitialized) {
		auto signalHandler = [](int signum ) {
			setRunning(false);
		};
		signal(SIGINT,signalHandler);
		handlerInitialized = true;
	}
}

} /* Aux */
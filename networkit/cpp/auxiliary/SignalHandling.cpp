#include "SignalHandling.h"
#include <exception>
#include <atomic>
#include <signal.h>

namespace Aux {

namespace SignalHandling {

namespace {
	std::atomic<bool> receivedSIGINT(false);
	std::atomic<bool> rootSet(false);
	std::atomic<uint64_t> root(0);
	void sigHandler(int sig) {
		switch (sig) {
			case SIGINT: receivedSIGINT = true;
			default: break;
		}
	};
}

bool gotSIGINT() {
	return receivedSIGINT;
}

void setSIGINT(bool received) {
	receivedSIGINT = received;
}

void init(uint64_t caller) {
	if (!rootSet) {
		root = caller;
		signal(SIGINT,sigHandler);
		rootSet = true;
	}
}

void reset(uint64_t caller) {
	if (root == caller) {
		rootSet = false;
		receivedSIGINT = false;
		root = 0;
		signal(SIGINT,0);
	}
}

} /* SignalHandling */

} /* Aux */
#include "SignalHandling.h"
#include <exception>
#include <atomic>
#include <signal.h>

namespace Aux {

namespace SignalHandling {

namespace {
	std::atomic<bool> receivedSIGINT(false);
	bool handlerInitialized = false;
	uint64_t root = 0;
	bool rootSet = false;
	void sigHandler(int sig) {
		switch (sig) {
			case SIGINT: receivedSIGINT.store(true);
			default: break;
		}
	};
	//struct sigaction act;
	//memset (&act, 0, sizeof (act));
	//act.sa_handler = sigHandler;
}

bool gotSIGINT() {
	return receivedSIGINT.load();
}

void setSIGINT(bool received) {
	receivedSIGINT.store(received);
}

void signalHandler(int signum) {
	setSIGINT(true);
}

void init(uint64_t caller) {
	if (!rootSet) {
		root = caller;
		if (!handlerInitialized) {
			//sigaction(SIGINT,&act,0);
			signal(SIGINT,sigHandler);
		}
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
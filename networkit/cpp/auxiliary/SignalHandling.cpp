#include "SignalHandling.h"
#include <exception>
#include <atomic>
#include <signal.h>

namespace Aux {

namespace SignalHandling {

namespace {
	std::atomic<bool> receivedSIGINT(false);
	std::atomic<bool> rootSet(false);
	std::atomic<SignalHandler*> root(0);
	std::atomic<void (*)(int)> prev_handler;
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

void init(SignalHandler* caller) {
	#pragma omp critical (SignalHandlerCritical)
	{
		if (!rootSet) {
			root = caller;
			prev_handler = signal(SIGINT,sigHandler);
			rootSet = true;
		}
	}
}

void reset(SignalHandler* caller) {
	#pragma omp critical (SignalHandlerCritical)
	{
		if (root == caller) {
			rootSet = false;
			receivedSIGINT = false;
			root = 0;
			signal(SIGINT,prev_handler);
		}
	}
}


} /* SignalHandling */

SignalHandler::SignalHandler() {
	Aux::SignalHandling::init(this);
}

SignalHandler::~SignalHandler() {
	Aux::SignalHandling::reset(this);
}

void SignalHandler::assureRunning() {
	if (Aux::SignalHandling::gotSIGINT()) {
		throw Aux::SignalHandling::InterruptException();
	}
}

bool SignalHandler::isRunning() {
	return !Aux::SignalHandling::gotSIGINT();
}

} /* Aux */
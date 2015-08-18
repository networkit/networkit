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
	if (!rootSet) {
		root = caller;
		signal(SIGINT,sigHandler);
		rootSet = true;
	}
}

void reset(SignalHandler* caller) {
	if (root == caller) {
		rootSet = false;
		receivedSIGINT = false;
		root = 0;
		signal(SIGINT,0);
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
#include <csignal>
#include <networkit/GlobalState.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>

namespace Aux {

namespace {

// this function is the signal handler
void sigHandler(int sig) {
    switch (sig) {
    case SIGINT:
        NetworKit::GlobalState::setReceivedSIGINT(true);
    default:
        break;
    }
}
} // namespace

SignalHandler::SignalHandler() {

#pragma omp critical(SignalHandlerCritical)
    {
        if (!NetworKit::GlobalState::getRootSet()) {
            // Registers the function `sigHandler` to be invoked, when CTRL+C/SIGINT is received.
            // Also registers `this` to be the owner of the signal handling resource (via RAII
            // responsible for cleanup and unregistering)
            NetworKit::GlobalState::setRoot(this);
            NetworKit::GlobalState::setPrevHandler(signal(SIGINT, sigHandler));
            NetworKit::GlobalState::setRootSet(true);
        }
    }
}

SignalHandler::~SignalHandler() {
    // Resets receivedSIGINT to false and rootSet to false to allow a new initialization,
    // if `this` is the root.
#pragma omp critical(SignalHandlerCritical)
    {
        if (NetworKit::GlobalState::getRoot() == this) {
            NetworKit::GlobalState::setRootSet(false);
            NetworKit::GlobalState::setReceivedSIGINT(false);
            NetworKit::GlobalState::setRoot(0);
            signal(SIGINT, NetworKit::GlobalState::getPrevHandler());
        }
    }
}

void SignalHandler::assureRunning() {
    if (NetworKit::GlobalState::getReceivedSIGINT()) {
        throw InterruptException();
    }
}

bool SignalHandler::isRunning() {
    return !NetworKit::GlobalState::getReceivedSIGINT();
}

} // namespace Aux

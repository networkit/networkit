#include <atomic>
#include <csignal>
#include <exception>
#include <networkit/GlobalState.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>

namespace Aux {

namespace SignalHandling {

namespace {
void sigHandler(int sig) {
    switch (sig) {
    case SIGINT:
        NetworKit::GlobalState::setReceivedSIGINT(true);
    default:
        break;
    }
}
} // namespace

void init(SignalHandler *caller) {
#pragma omp critical(SignalHandlerCritical)
    {
        if (!NetworKit::GlobalState::getRootSet()) {
            NetworKit::GlobalState::setRootSet(caller);
            NetworKit::GlobalState::setPrevHandler(signal(SIGINT, sigHandler));
            NetworKit::GlobalState::setRootSet(true);
        }
    }
}

void reset(SignalHandler *caller) {
#pragma omp critical(SignalHandlerCritical)
    {
        if (NetworKit::GlobalState::getRoot() == caller) {
            NetworKit::GlobalState::setRootSet(false);
            NetworKit::GlobalState::setReceivedSIGINT(false);
            NetworKit::GlobalState::setRoot(0);
            signal(SIGINT, NetworKit::GlobalState::getPrevHandler());
        }
    }
}

} // namespace SignalHandling

SignalHandler::SignalHandler() {
    Aux::SignalHandling::init(this);
}

SignalHandler::~SignalHandler() {
    Aux::SignalHandling::reset(this);
}

void SignalHandler::assureRunning() {
    if (NetworKit::GlobalState::getReceivedSIGINT()) {
        throw Aux::SignalHandling::InterruptException();
    }
}

bool SignalHandler::isRunning() {
    return !NetworKit::GlobalState::getReceivedSIGINT();
}

} // namespace Aux

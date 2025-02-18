#ifndef NETWORKIT_AUXILIARY_SIGNAL_HANDLING_HPP_
#define NETWORKIT_AUXILIARY_SIGNAL_HANDLING_HPP_

#include <cstdint>
#include <exception>

namespace Aux {

/**
 * Signal handler class used to support CTRL+c user interrupts.
 * This class works by registering a interrupt handler during its lifetime; to check if an interrupt
 * has been received or to immediately stop running a (long-running) algorithm, use the
 * `assureRunning` or `isRunning` member functions.
 * @note Make sure to delete this object when it is no longer required to free up the interrupt
 * handler again for other code such as the python runtime!
 */
class SignalHandler {
public:
    SignalHandler();

    ~SignalHandler();

    void assureRunning();

    bool isRunning();

    /**
     * Special Exception to indicate that a SIGINT has been received.
     */
    class InterruptException : public std::exception {
    public:
        InterruptException() : std::exception() {}
        const char *what() const noexcept override { return "Received CTRL+C/SIGINT"; }
    };
};

} // namespace Aux

#endif // NETWORKIT_AUXILIARY_SIGNAL_HANDLING_HPP_

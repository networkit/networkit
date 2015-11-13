#ifndef SIGNALHANDLING_H
#define SIGNALHANDLING_H

#include <cstdint>
#include <exception>

namespace Aux {

	class SignalHandler {
	public:
		SignalHandler();

		~SignalHandler();

		void assureRunning();

		bool isRunning();
	};
	
namespace SignalHandling {
		/**
		 * Special Exception to indicate, that a SIGINT has been received.
		 */
		class InterruptException : public std::exception {
		public:
			InterruptException() : std::exception() {}
			virtual const char* what() const noexcept
			{
			    return "Received CTRL+C/SIGINT";
			}
		};

		/**
		 * Returns true, if CTRL+C/SIGINT has been received, false otherwise.
		 * @return A boolean indicating if SIGINT has been received since construction/reset
		 */
		bool gotSIGINT();
		
		/**
		 * Sets the value of `receivedSIGINT`
		 * @param received The new value for `receivedSIGINT`
		 */
		void setSIGINT(bool received);

		/**
		 * Registers the function `setRunning` to be invoked, when CTRL+C/SIGINT is received.
		 * Also registers `caller` to be the root calling object, if none is registered yet.
		 * @param caller address of the calling object
		 */
		void init(SignalHandler* caller);

		/**
		 * Resets receivedSIGINT to false and rootSet to false to allow a new initialization, 
		 * if the calling object is the same as the root.
		 * @param caller address of the calling object
		 */
		void reset(SignalHandler* caller);
	};

}

#endif /* SIGNALHANDLING_H */
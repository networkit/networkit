#ifndef SIGNALHANDLING_H
#define SIGNALHANDLING_H

#include <csignal>
#include <cstdint>
#include <exception>
//#include <atomic>

namespace Aux {
	/**
	 * Special Exception to indicate, that a SIGINT has been received.
	 */
	class InterruptException : public std::exception {
	public:
		virtual const char* what() const throw()
		{
		    return "Received CTRL+C/SIGINT";
  		}
	};

	class SignalHandling {
	private:
		//static std::atomic<bool> running;
		//static std::atomic<bool> handlerInitialized;
		static bool running;
		static bool handlerInitialized;
		static bool rootSet;
		static uint64_t root;

	public:
		/**
		 * Returns true, as long as no CTRL+C/SIGINT has been received, 
		 * throws an InterruptExcpetion otherwise.
		 * @return A boolean indicating whether an algorithm still may run.
		 */
		static bool isRunning();
		
		/**
		 * Sets the value of `running`
		 * @param running The new value for `running`
		 */
		static void setRunning(bool running);

		/**
		 * Registers the function `setRunning` to be invoked, when CTRL+C/SIGINT is received.
		 * Also registers `caller` to be the root calling object, if none is registered yet.
		 * @param caller address of the calling object
		 */
		static void init(uint64_t caller);

		/**
		 * Resets running to true and rootSet to false to allow a new initialization, 
		 * if the calling object is the same as the root.
		 * @param caller address of the calling object
		 */
		static void reset(uint64_t caller);
	};
}

#endif /* SIGNALHANDLING_H */
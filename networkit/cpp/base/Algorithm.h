#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <string>
#include <stdexcept>

namespace NetworKit {
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

class Algorithm {
protected:
	/**
	 * The method that actually holds the implementation of an algorithm. Subclasses need to override it.
	 */
	virtual void runImpl() = 0;
	/**
	 * A boolean variable indicating whether an algorithm has finished its computation or not.
	 */
	bool hasRun;
public:
	/**
	 * Constructor to the algorithm base class.
	 */
	Algorithm();

	/**
	 * Virtual default destructor
	 */
	virtual ~Algorithm() = default;

	/**
	 * The generic run method which calls runImpl() and takes care of setting @hasRun to the appropriate value.
	 */
	void run();

	/**
	 * Uses SignalHandling to check whether a SIGINT signal has been received. If so, it throws an InterruptException.
	 */
	void assureRunning();

	/**
	 * Indicates whether an algorithm has completed computation or not.
	 * @return The Value of @hasRun.
	 */
	bool hasFinished() const;

	/**
	 * Returns a string with the algorithm's name and its parameters, if there are any. Subclasses should override it.
	 * @return The string representation of the algorithm.
	 */
	virtual std::string toString() const;
};

} /* NetworKit */

#endif /* ALGORITHM_H */
#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <string>
#include <stdexcept>

namespace NetworKit {
	/**
	 * Special Exception to indicate, that a SIGINT has been received.
	 */
	class InterruptException : public std::runtime_error {
	public:
		InterruptException() : std::runtime_error("") {}
		virtual const char* what() const noexcept
		{
		    return "Received CTRL+C/SIGINT";
  		}
	};

class Algorithm {
protected:
	virtual void runImpl() = 0;
	bool hasRun;
public:
	Algorithm();
	void run();
	void assureRunning();
	bool hasFinished() const;
	virtual std::string toString() const;
};

} /* NetworKit */

#endif /* ALGORITHM_H */
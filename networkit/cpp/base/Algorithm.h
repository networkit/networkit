#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <string>
#include <stdexcept>

namespace NetworKit {

class Algorithm {
protected:
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
	virtual void run() = 0;

	/**
	 * Indicates whether an algorithm has completed computation or not.
	 * @return The Value of @hasRun.
	 */
	bool hasFinished() const {
		return hasRun;
	};

	/**
	 * Assure that the algorithm has been run, throws a std::runtime_error otherwise.
	 */
	void assureFinished() const {
		if (!hasRun) throw std::runtime_error("Error, run must be called first");
	};

	/**
	 * Returns a string with the algorithm's name and its parameters, if there are any. Subclasses should override it.
	 * @return The string representation of the algorithm.
	 */
	virtual std::string toString() const;

	/**
	 * @return True if algorithm can run multi-threaded. 
	 */
	virtual bool isParallel() const;
};

} /* NetworKit */

#endif /* ALGORITHM_H */

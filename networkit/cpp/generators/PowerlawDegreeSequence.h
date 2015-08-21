#ifndef POWERLAWDEGREESEQUENCE_H
#define POWERLAWDEGREESEQUENCE_H

#include "../Globals.h"
#include <vector>

namespace NetworKit {

class PowerlawDegreeSequence {
public:
	/**
	 * Generates a powerlaw degree sequence with the given minimum and maximum degree, the powerlaw exponent gamma.
	 *
	 * @param minDeg   The minium degree
	 * @param maxDeg   The maximum degree
	 * @param gamma    The powerlaw exponent
	 */
	PowerlawDegreeSequence(count minDeg, count maxDeg, double gamma);

	/**
	 * Tries to set the minimum degree such that the specified average degree is expected.
	 *
	 * @throws std::runtime_error If it is not possible to find a minimum degree such that the expected average is @a avgDeg.
	 * @param avgDeg   The average degree
	 */
	void setMinimumFromAverageDegree(double avgDeg);

	/**
	 * Sets the minimum degree.
	 *
	 * @param minDeg The degree that shall be set as minimum degree
	 */
	void setMinimumDegree(count minDeg);

	/**
	 * Gets the minimum degree.
	 *
	 * @return The minimum degree.
	 */
	count getMinimumDegree() const;

	/**
	 * Execute the generation process
	 */
	void run();

	/**
	 * Returns the expected average degree.
	 *
	 * @return The expected average degree.
	 */
	double getExpectedAverageDegree() const;

	/**
	 * Returns a degree sequence of @a numNodes degrees with even degree sum.
	 *
	 * @param numNodes The number of nodes
	 * @return The generated degree sequence.
	 */
	std::vector<count> getDegreeSequence(count numNodes) const;

	/**
	 * Returns a degree drawn at random with a power law distribution
	 *
	 * @return A degree that follows the generated distribution.
	 */
	count getDegree() const;
private:
	count minDeg, maxDeg;
	double gamma;
	bool hasRun;
	std::vector<double> cumulativeProbability;
};

} // namespace NetworKit

#endif // POWERLAWDEGREESEQUENCE_H

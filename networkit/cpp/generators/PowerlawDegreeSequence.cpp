/*
 *
 */

#include <cmath>
#include <algorithm>

#include "PowerlawDegreeSequence.h"
#include "../auxiliary/Random.h"
#include "../auxiliary/NumericTools.h"

NetworKit::PowerlawDegreeSequence::PowerlawDegreeSequence(NetworKit::count minDeg, NetworKit::count maxDeg, double gamma) :
	minDeg(minDeg), maxDeg(maxDeg), gamma(gamma), hasRun(false) {
	if (minDeg > maxDeg) throw std::runtime_error("Error: minDeg must not be larger than maxDeg");
	if (gamma > -1) throw std::runtime_error("Error: gamma must be lower than -1");
}

void NetworKit::PowerlawDegreeSequence::setMinimumDegree(NetworKit::count minDeg) {
	this->minDeg = minDeg;
	hasRun = false;
}

void NetworKit::PowerlawDegreeSequence::setMinimumFromAverageDegree(double avgDeg) {
	count dmin_l = 1;
	count dmin_r = maxDeg;
	setMinimumDegree(dmin_l); run();
	double average_l = getExpectedAverageDegree();
	double average_r = maxDeg;

	if (average_l > avgDeg) {
		throw std::runtime_error("The average degree is too low");
	}

	if (average_r < avgDeg) {
		throw std::runtime_error("The average degree must not be higher than the maximum degree");
	}

	while (dmin_l + 1 < dmin_r) {
		setMinimumDegree((dmin_r + dmin_l) * 0.5); run();
		double avg = getExpectedAverageDegree();

		TRACE("Trying minDeg ", minDeg, ", this gives average ", avg, ", which should be between ", average_r, " and ", average_l);

		if (avg > avgDeg) {
			average_r = avg;
			dmin_r = minDeg;
		} else {
			average_l = avg;
			dmin_l = minDeg;
		}
	}

	if (avgDeg - average_l < average_r - avgDeg) {
		minDeg = dmin_l;
	} else {
		minDeg = dmin_r;
	}

	hasRun = false;
}

NetworKit::count NetworKit::PowerlawDegreeSequence::getMinimumDegree() const {
	return minDeg;
}

void NetworKit::PowerlawDegreeSequence::run() {
	cumulativeProbability.clear();
	cumulativeProbability.reserve(maxDeg - minDeg + 1);

	double sum = 0;

	for (double d = minDeg; d <= maxDeg; ++d) {
		sum += std::pow(d, gamma);
		cumulativeProbability.push_back(sum);
	}

	for (double & prob : cumulativeProbability) {
		prob /= sum;
	}

	hasRun= true;
}

double NetworKit::PowerlawDegreeSequence::getExpectedAverageDegree() const {
	if (!hasRun) throw std::runtime_error("Error: run needs to be called first");

	double average = cumulativeProbability[0] * minDeg;
	for (count i = 1; i < cumulativeProbability.size(); ++i) {
		average += (cumulativeProbability[i] - cumulativeProbability[i-1]) * (i + minDeg);
	}

	return average;
}

std::vector< NetworKit::count > NetworKit::PowerlawDegreeSequence::getDegreeSequence(NetworKit::count numNodes) const {
	std::vector<count> degreeSequence;

	if (!hasRun) throw std::runtime_error("Error: run needs to be called first");

	degreeSequence.reserve(numNodes);
	count degreeSum = 0;

	for (count i = 0; i < numNodes; ++i) {
		degreeSequence.push_back(getDegree());
		degreeSum += degreeSequence.back();
	}

	if (degreeSum % 2 != 0) {
		(*std::max_element(degreeSequence.begin(), degreeSequence.end()))--;
	}

	return degreeSequence;
}

NetworKit::count NetworKit::PowerlawDegreeSequence::getDegree() const {
	if (!hasRun) throw std::runtime_error("Error: run needs to be called first");

	return std::distance(cumulativeProbability.begin(), std::lower_bound(cumulativeProbability.begin(), cumulativeProbability.end(), Aux::Random::probability())) + minDeg;
}

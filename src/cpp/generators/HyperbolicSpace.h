/*
 * HyperbolicSpace.h
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#ifndef HYPERBOLICSPACE_H_
#define HYPERBOLICSPACE_H_

namespace NetworKit {

class HyperbolicSpace {
public:
	HyperbolicSpace();
	virtual ~HyperbolicSpace();
	HyperbolicSpace(double R);
	static double getDistance(double firstangle, double firstR, double secondangle, double secondR);//TODO: replace with coordinates
	static double getDistancePrecached(double firstangle, double firstRcosh, double firstRsinh, double secondangle, double secondRcosh, double secondRsinh);
	double getRadius();


private:
	double radius;
	//one could add some caching here, but this should be done properly in an object. Well, or maybe not.
	static double lastR;
	static double coshlastR;
	static double sinhlastR;
};
}

#endif /* HYPERBOLICSPACE_H_ */

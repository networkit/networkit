/*
 * DynamicHyperbolicGenerator.h
 *
 *  Created on: 29.07.2014
 *      Author: moritzl
 */

#ifndef DYNAMICHYPERBOLICGENERATOR_H_
#define DYNAMICHYPERBOLICGENERATOR_H_

namespace NetworKit {

class DynamicHyperbolicGenerator {
public:
	DynamicHyperbolicGenerator();
	virtual ~DynamicHyperbolicGenerator();

private:
	double factorgrowth;
	double moveDistance;
	double spread;
	double initialfactor;

};

} /* namespace NetworKit */
#endif /* DYNAMICHYPERBOLICGENERATOR_H_ */

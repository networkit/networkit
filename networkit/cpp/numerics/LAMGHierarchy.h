/*
 * LAMGHierarchy.h
 *
 *  Created on: 14.11.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef LAMGHIERARCHY_H_
#define LAMGHIERARCHY_H_

#include "MultigridHierarchy.h"

namespace NetworKit {

class LAMGHierarchy : public MultigridHierarchy {
public:
	enum Type{ELIMINATION, AGGREGATION};

	LAMGHierarchy() {}

	void addLevel(const Matrix &laplacian, const Matrix &interpolationMatrix, Type type);
	using MultigridHierarchy::addLevel;

	Type getType(const index level) const;

private:
	std::vector<Type> type;
};

} /* namespace NetworKit */

#endif /* LAMGHIERARCHY_H_ */

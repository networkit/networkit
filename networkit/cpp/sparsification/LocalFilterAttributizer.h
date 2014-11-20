/*
 * LocalFilterAttributizer.h
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann, Gerd Lindner
 */

#ifndef LOCALLOGATTRIBUTIZER_H
#define LOCALLOGATTRIBUTIZER_H

#include "../edgeattributes/EdgeAttribute.h"

namespace NetworKit {

template<typename InType>
class LocalFilterAttributizer : public EdgeAttribute<double> {

public:
	LocalFilterAttributizer(const Graph &graph, const std::vector< InType > &attribute, bool logarithmic = true, bool bothRequired = false);
	virtual std::vector< double > getAttribute() override;

private:
	const Graph& graph;
	const std::vector<InType>& attribute;
	bool bothRequired;
	bool logarithmic;

};
} // namespace NetworKit

#endif // LOCALLOGATTRIBUTIZER_H

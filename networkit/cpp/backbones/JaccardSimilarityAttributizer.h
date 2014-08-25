#ifndef JACCARDSIMILARITYATTRIBUTIZER_H
#define JACCARDSIMILARITYATTRIBUTIZER_H

#include "AttributeGenerator.h"

namespace NetworKit {

class JaccardSimilarityAttributizer : public AttributeGenerator<int, double> {
public:
	virtual std::vector<double> getAttribute(const Graph& g, const std::vector<int>& triangles) override;
};

}

#endif // JACCARDSIMILARITYATTRIBUTIZER_H

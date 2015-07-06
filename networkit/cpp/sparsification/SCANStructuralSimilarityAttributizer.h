
#ifndef SCANSTRUCTURALSIMILARITY_H
#define SCANSTRUCTURALSIMILARITY_H

#include "../edgeattributes/EdgeAttribute.h"

namespace NetworKit {

	class SCANStructuralSimilarityAttributizer : public EdgeAttribute<double> {

	public:
		SCANStructuralSimilarityAttributizer(const Graph& graph, const std::vector<count>& triangles);
		virtual std::vector<double> getAttribute() override;

	private:
		const Graph& graph;
		const std::vector<count>& triangles;

	};

} // namespace NetworKit

#endif // SCANSTRUCTURALSIMILARITY_H

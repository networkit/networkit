
#ifndef SCANSTRUCTURALSIMILARITY_H
#define SCANSTRUCTURALSIMILARITY_H

#include "../edgescores/EdgeScore.h"

namespace NetworKit {

	class SCANStructuralSimilarityScore : public EdgeScore<double> {

	public:
		SCANStructuralSimilarityScore(const Graph& G, const std::vector<count>& triangles);
		virtual void run() override;
		virtual double score(edgeid eid) override;
		virtual double score(node u, node v) override;

	private:
		const std::vector<count>& triangles;

	};

} // namespace NetworKit

#endif // SCANSTRUCTURALSIMILARITY_H

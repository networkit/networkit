#ifndef COVERF1SIMILARITY_H_
#define COVERF1SIMILARITY_H_

#include "LocalCoverEvaluation.h"

namespace NetworKit {

/**
 * @ingroup community
 * Compare a given cover to a reference cover using the F1 measure.
 * This is a typical similarity measure used to compare the found
 * overlapping community structure to a ground truth community
 * structure. Each cluster is compared to the best-matching reference
 * cluster (in terms of highest F1 score). A value of 1 indicates
 * perfect agreement while a while of 0 indicates complete
 * disagreement. An example where this measure is used is the
 * following paper:
 *
 * Alessandro Epasto, Silvio Lattanzi, and Renato Paes
 * Leme. 2017. Ego-Splitting Framework: from Non-Overlapping to
 * Overlapping Clusters. In Proceedings of the 23rd ACM SIGKDD
 * International Conference on Knowledge Discovery and Data Mining
 * (KDD '17). ACM, New York, NY, USA, 145-154. DOI:
 * https://doi.org/10.1145/3097983.3098054
 */
class CoverF1Similarity : public LocalCoverEvaluation {
public:
	/**
	 * Initialize the cover F1 similarity.
	 *
	 * @param G The graph on which the evaluation is performed.
	 * @param C The cover that shall be evaluated.
	 * @param reference The reference cover to which @a C shall be compared.
	 */
	CoverF1Similarity(const Graph& G, const Cover& C, const Cover& reference);

	/**
	 * Execute the algorithm.
	 */
	virtual void run() override;

	/**
	 * @return false - smaller is not better, larger values indicate better matching clusters.
	 */
	virtual bool isSmallBetter() const override { return false; }

	/**
	 * @return false - this algorithm has not been parallelized.
	 */
	virtual bool isParallel() const override { return false; }
private:
	const Cover &reference;
};

}

#endif

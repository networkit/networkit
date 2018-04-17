#ifndef TIM_H_
#define TIM_H_

#include "../base/Algorithm.h"
#include "../graph/Graph.h"

#include <unordered_set>

namespace NetworKit {

/**
 * @ingroup influence
 * Two-phase influence maximization algorithm by Tang et al.
 */
class TwoPhaseInfluenceMaximization final : public Algorithm {

public:

    enum class Model {
        INDEPENDENT_CASCADE,
        LINEAR_THRESHOLD
    };


    /**
     * Constructs an instance of the two-phase influence maximization algorithm.
     *
     * This instance will find an (1 - 1/e - @a epsilon)-approximate solution to the influence maximization problem with
     * initial seed set size @a k on graph @a G with a probability of at least (1 - 2 * n ^ - @a l).
     *
     * The expected runtime of the algorithm is O((@a k + @a l) (m + n) log n / (@a epsilon ^ 2)).
     */
    TwoPhaseInfluenceMaximization(const Graph& G, count k, double epsilon, double l,
            Model model = Model::INDEPENDENT_CASCADE);

    void run() override;

    std::string toString() const override {
        return "Two-phase influence maximization algorithm by Tang et al.";
    }

    bool isParallel() const override {
        return false;
    }

    std::unordered_set<node> topKInfluencers() {
        return influencers;
    }

private:

    const Graph& G;
    const count k;
    const double epsilon;
    const double l;
    const Model model;

    std::unordered_set<node> influencers{};

};

}

#endif /* TIM_H_ */
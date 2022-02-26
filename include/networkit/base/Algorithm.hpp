
#ifndef NETWORKIT_BASE_ALGORITHM_HPP_
#define NETWORKIT_BASE_ALGORITHM_HPP_

#include <stdexcept>
#include <string>

namespace NetworKit {

class Algorithm {
protected:
    bool hasRun = false;

public:
    virtual ~Algorithm() = default;

    /**
     * The generic run method which calls runImpl() and takes care of setting @ref hasRun to the
     * appropriate value.
     */
    virtual void run() = 0;

    /**
     * Indicates whether an algorithm has completed computation or not.
     *
     * @return The value of @ref hasRun.
     */
    bool hasFinished() const noexcept { return hasRun; }

    /**
     * Assure that the algorithm has been run, throws a std::runtime_error otherwise.
     */
    void assureFinished() const {
        if (!hasRun)
            throw std::runtime_error("Error, run must be called first");
    }
};

} // namespace NetworKit

#endif // NETWORKIT_BASE_ALGORITHM_HPP_

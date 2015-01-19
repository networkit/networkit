#ifndef ALGORITHM_H
#define ALGORITHM_H

namespace NetworKit {
class Algorithm {
protected:
	virtual void runImpl() = 0;
	bool hasRun;
public:
	Algorithm();
	void run();
};

} /* NetworKit */

#endif /* ALGORITHM_H */
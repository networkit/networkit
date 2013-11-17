/*
 * ShellList.h
 *
 *  Created on: Nov 4, 2013
 *      Author: Lukas Barth
 */
#ifndef SHELLLIST_H_
#define SHELLLIST_H_

#include <list>
#include <vector>
#include <map>
#include "../graph/Graph.h"

using namespace NetworKit;

namespace Aux {

/**
 * Noise is random addition to a signal. This class provides methods
 * which add random numbers to their inputs in order to enable randomization.
 */
class ShellList {

private:
	const Graph *g;

	std::vector<std::list<node>> shells;
	std::map<node, count> nodeToShell;
	std::map<node, std::list<node>::iterator> listHandle;

protected:

public:
	ShellList(const Graph* orig_g);

	void insert(node n);
	void decreaseDegree(node n);
	count getCurrentShell(node n);
	std::list<node>::iterator getShelliterator(count shell);
	count size();

	void forEachNodeInShell(count shell, std::function<void (node)> func);

	std::vector<count> getCoreness();
};

}
#endif
/*
 * ShellList.h
 *
 * FIXME: this does not belong in the Aux namespace
 *
 *  Created on: Nov 4, 2013
 *      Author: Lukas Barth, David Weiß
 */
#ifndef SHELLLIST_H_
#define SHELLLIST_H_

#include <list> // std::list, doubly linked list
#include <vector> // std::vector, dynamic array
#include <functional> // std::function, typesafe function pointer
#include "../graph/Graph.h"

using namespace NetworKit;

namespace Aux {

class ShellList {

public:
	ShellList(const std::vector<count>& seq);
  
  /**
   * Moves the node one shell downwards.
   */
	void decreaseShell(const node v);
  
  /**
   * Executes func for each node in the given shell.
   */
	void forEachNodeInShell(const count shell, std::function<void (node)> func);
  
  /**
   * @return Shell the given node is currently assigned to. After initialization this equals its degree.
   */
	count getShell(const node v) const
  {
    return this->nodeShell[v];
  }
  
  /**
   * @return Vector with the shells of all nodes indexed by node.
   */
	std::vector<count> getShells() const
  {
    return this->nodeShell;
  }


  /**
   * @return Size of the shell list.
   * Note that we don't clean up empty shells, so this will always be equal
   * to the maximum degree of the graph on initialization.
   */
	count size() const
  {
    return this->shells.size();
  }

	bool isShellEmpty(count shell) const {
		return shells[shell].empty();
	}

	node popVertexOfShell(count shell) {
		assert(! shells[shell].empty());
		node result = (* shells[shell].begin());
		shells[shell].pop_front();
		return result;
	}


private:
	std::vector<count> seq;

	std::vector<std::list<node>> shells;
	std::vector<std::list<node>::iterator> nodeHandle;
	std::vector<count> nodeShell;
};

}//namespace Aux

#endif//SHELLLIST_H_

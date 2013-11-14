/*
 * ShellList.cpp
 *
 *  Created on: Nov 4, 2013
 *      Author: Lukas Barth, David Weiß
 */

#include <algorithm> // std::for_each, loop iterating over containers and executing a given function for each element
#include "ShellList.h"

namespace Aux {

ShellList::ShellList(const Graph *graph):
	g(graph),
	shells(std::vector<std::list<node>>(graph->maxDegree(), std::list<node>())),
  nodeHandle(std::vector<count>(graph->numberOfNodes())),
  nodeShell(std::vector<count>(graph->numberOfNodes()))
{
 	this->g->forNodes([&](node v) {
    int shell = this->g->degree(v)
    this->shells[shell].push_front(v);
    this->nodeHandle[v] = this->shells[shell].begin();
    this->nodeShell[v] = shell;
 	});
}

void 
ShellList::decreaseDegree(const node v) 
{
  count shell = this->nodeShell[v];
  assert(shell > 0);
  
	this->shells[shell].erase(this->nodeHandle[v]);
  
  --this->nodeShell[v];
  --shell;
  
	this->shells[shell].push_back(v);
	this->nodeHandle[v] = this->shells[shell].end();
	--this->nodeHandle[v]; // iterator arithmetic, necessary because .end() does not point on the last element but after it
}

void 
ShellList::forEachNodeInShell(const count shell, std::function<void (node)> func)
{
  assert(shell >= 0);
  assert(shell < this->shells.size());
	std::for_each(this->shells[shell].begin(), this->shells[shell].end(), func);
}

}//namespace Aux
/*
 * ShellList.cpp
 *
 *  Created on: Nov 4, 2013
 *      Author: Lukas Barth, David Weiß
 */

#include <algorithm> // std::for_each, loop iterating over containers and executing a given function for each element
#include "ShellList.h"
#include "../auxiliary/Log.h"

namespace Aux {

ShellList::ShellList(const std::vector<count>& sequence): seq(sequence)
{
	assert(! seq.empty());
	count n = seq.size();
	count k = (* std::max_element(seq.begin(), seq.end())) + 1;
	shells.resize(k);
	nodeHandle.resize(n);
	nodeShell.resize(n);

	DEBUG("allocation in shell list finished");

	for (uint64_t i = 0; i < seq.size(); ++i) {
 		count shell = seq[i];
 		this->shells[shell].push_front(i);
 		this->nodeHandle[i] = this->shells[shell].begin();
 		this->nodeShell[i] = shell;
	}
}

void 
ShellList::decreaseShell(const node v) 
{
	count shell = this->nodeShell[v];
	assert(shell > 0);
  
	this->shells[shell].erase(this->nodeHandle[v]);
  
	--this->nodeShell[v];
	--shell;
  
	this->shells[shell].push_back(v);
	this->nodeHandle[v] = std::prev(this->shells[shell].end());
}

void 
ShellList::forEachNodeInShell(const count shell, std::function<void (node)> func)
{
	assert(shell >= 0);
	assert(shell < this->shells.size());
	std::for_each(this->shells[shell].begin(), this->shells[shell].end(), func);
}

}//namespace Aux

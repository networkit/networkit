/*
 * ShellList.cpp
 *
 *  Created on: Nov 4, 2013
 *      Author: Lukas Barth
 */

#include "ShellList.h"

namespace Aux {

ShellList::ShellList(const Graph *orig_g):
	g(orig_g),
	shells(std::vector<std::list<node>>(orig_g->numberOfNodes(), std::list<node>()))
	{
}

void
ShellList::insert(node n) 
{
	this->shells[this->g->degree(n)].push_front(n);
	this->nodeToShell[n] = this->g->degree(n);
	this->listHandle[n] = this->shells[this->g->degree(n)].begin();
}

void 
ShellList::decreaseDegree(node n) 
{
	std::list<node>::iterator nIt = this->listHandle[n];
	count curDeg = this->nodeToShell[n];

	this->shells[curDeg].erase(nIt);
	this->shells[curDeg-1].push_back(n);

	this->nodeToShell[n] -= 1;
	this->listHandle[n] = this->shells[curDeg-1].end();
	--this->listHandle[n];
}

count
ShellList::getCurrentShell(node n) 
{
	return this->nodeToShell[n];
}

std::list<node>::iterator
ShellList::getShelliterator(count shell) 
{
	return this->shells[shell].begin();
}

count
ShellList::size()
{
	return this->shells.size();
}

void 
ShellList::forEachNodeInShell(count shell, std::function<void (node)> func)
{
	for_each(this->shells[shell].begin(), this->shells[shell].end(), func);
}

std::vector<count>
ShellList::getCoreness()
{
	std::vector<count> coreness(this->nodeToShell.size());

	for (auto it = this->nodeToShell.begin(); it != this->nodeToShell.end(); ++it) {
		coreness[it->first] = it->second;
	}

	return coreness;
}

}
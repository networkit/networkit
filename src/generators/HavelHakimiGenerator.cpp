/*
 * HavelHakimiGenerator.cpp
 *
 *  Created on: Dec 10, 2013
 *      Author: Henning
 *      Contributors: Hoske/Weisbarth
 */

#include "HavelHakimiGenerator.h"

namespace NetworKit {

HavelHakimiGenerator::HavelHakimiGenerator(const std::vector<count>& sequence) :
		seq(sequence) {
	std::sort(seq.begin(), seq.end(), std::greater<count>());
	realizable = isRealizable();
}

bool HavelHakimiGenerator::isRealizable() {
	DEBUG("check if sequence is realizable");
	realizable = true;
	count n = seq.size();

	/* First inequality. */
	count deg_sum = 0;
	for (count i = 0; i < n; ++i) {
		if (seq[i] < 0 || seq[i] >= n) {
			realizable = false;
			DEBUG("not realizable: ", seq[i], ", n: ", n);
			break;
		}
		deg_sum += seq[i];
	}

	if (deg_sum % 2 != 0) {
		DEBUG("not realizable");
		realizable = false;
	}

	/* Second inequality. */
	deg_sum = 0;
	for (count j = 0; j < n; ++j) {
		deg_sum += seq[j];

		/* sum of min(deg(i), j) for i from j + 1 to n - 1. */
		count min_deg_sum = 0;
		for (count i = j + 1; i < n; ++i) {
			min_deg_sum += std::min(seq[i], j + 1);
		}

		if (deg_sum > (j + 1) * j + min_deg_sum) {
			DEBUG("not realizable");
			realizable = false;
			break;
		}
	}

	return realizable;
}

bool HavelHakimiGenerator::getRealizable() const {
	return realizable;
}

Graph HavelHakimiGenerator::generate() {
	count n = seq.size();
	count vol = std::accumulate(seq.begin(), seq.end(), 0);

	if (!realizable) {
		WARN("Degree sequence not realizable or not checked for realizability yet! Will return empty graph!");
		Graph G;
		return G;
	} else {
		DEBUG("Degree sequence is realizable, continue with generation algorithm.");
		Graph G(n);
		count numDegVals = (* std::max_element(seq.begin(), seq.end())) + 1;


		typedef std::pair<count, node> DeficitAndNode;
		typedef std::vector<std::list<DeficitAndNode> > Buckets;

		// put nodes in appropriate lists
		Buckets nodesByDeficit(n);
		std::list<DeficitAndNode>::iterator listIter;
		std::vector<std::list<DeficitAndNode>::iterator > nodePointer;

		for(node v = 0; v < n; v++) {
			nodesByDeficit[seq[v]].push_front(std::make_pair(seq[v], v));
			nodePointer.push_back(nodesByDeficit[seq[v]].begin());
		}

		index maxDeficit = numDegVals - 1;

		while (maxDeficit > 0) {
			// process node in largest bucket
			while(! nodesByDeficit[maxDeficit].empty()) {
				// get element
				listIter = nodesByDeficit[maxDeficit].begin();
				count deficit = listIter->first;
				node currentVertex = listIter->second;

				// delete it
				nodesByDeficit[maxDeficit].pop_front();
				std::set<node> toChange;

				// connect corresponding vertex with the following ones
				index currentNeighborList = maxDeficit;
				while (deficit > 0) {
					// search for candidates in current list
					for (auto elem : nodesByDeficit[currentNeighborList]) {
						// connect
						node nextNeighbor = elem.second;
						G.addEdge(currentVertex, nextNeighbor);
						--deficit;

						// mark as "to change"
						toChange.insert(nextNeighbor);
					}
					--currentNeighborList;
				}

				// remove entry for v, move connected vertices to appropriate lists
				for (auto v : toChange) {
					// delete from current list, insert to the lower one
					std::list<DeficitAndNode>::iterator iter = nodePointer[v];
					count neighDef = iter->first;
					nodesByDeficit[neighDef].erase(iter);
					--neighDef;
					if (neighDef > 0) {
						// insert adapted values
						iter->first = neighDef;
						nodesByDeficit[neighDef].push_front(* iter);
					}
				}
			}
			maxDeficit--;
		}

		//
		//
		////		/* Pairs of (remaining degree, node id). */
		////		std::vector<deficitAndNode> rem_deg(n);
		////		for (node v = 0; v < n; ++v) {
		////			rem_deg[v] = {seq[v], v};
		////		}
		//
		//		/* In each loop: fulfill degree requirement for one node. */
		//		for (count k = 0; k < n; ++k) {
		//			/* Distribute edges for node with largest remaining degree. */
		//
		//			// TODO: avoid explicit sorting by FM bucket data structure!!!
		//			std::sort(rem_deg.begin() + k, rem_deg.end(),
		//					std::greater<DeficitAndNode>());
		//			count& deg = rem_deg[k].first;
		//			node& v = rem_deg[k].second;
		//			assert(deg <= n - k - 1);
		//
		//			for (count l = k + 1; l <= k + deg; ++l) {
		//				rem_deg[l].first--;
		//				G.addEdge(v, rem_deg[l].second);
		//			}
//		}


//		std::vector<count> bucket(numDegVals, 0);
//		std::vector<count> sorted(n);
//
//		// count bucket entries
//		for (index i = 0; i < n; ++i) {
//			bucket[seq[i]]++;
//		}
//
//		// accumulate for descending order
//		count acc = 0;
//		for (index j = numDegVals-1; j >= 0; --j) {
//			acc += bucket[j];
//			bucket[j] = -bucket[j] + acc;
//		}
//
//		// sort
//		for (index i = 0; i < n; ++i) {
//			sorted[bucket[seq[i]]++] = i;
//		}
//
//		count current = 0; // vertex with largest deficit
//		count maxDeficit = numDegVals - 1; // largest active deficit
//
//		for (index e = vol / 2; e > 0; ) {
//			// pick the vertex with largest deficit
//			node v = sorted[current];
//
//			// its bucket number is its deficit
//
//			// add edges to the following vertices (deficit times)
//
//			// move vertices to appropriate buckets
//		}

#if 0
		 p = len(deg_sequence) ;
		 G=nx.empty_graph(p,create_using);
		 num_degs = [];
		 for currentDeficit in range(p):
				 num_degs.append([]);
		 dmax, dsum, n = 0, 0, 0;
		 for d in deg_sequence:
		 	 # Process only the non-zero integers
		 	 if d>0:
		 	 	 num_degs[d].append(n)
		 	 	 dmax, dsum, n = max(dmax,d), dsum+d, n+1

		 # Return graph if no edges
		 if n==0: return G;

		 modstubs = [(0,0)]*(dmax+1)
			# Successively reduce degree sequence by removing the maximum degree
		while n > 0:
			# Retrieve the maximum degree in the sequence
			while len(num_degs[dmax]) == 0:
				dmax -= 1;

		 	 	# Remove largest stub in list
		 	 	source = num_degs[dmax].pop()
		 	 	n -= 1

		 	 	# Reduce the next dmax largest stubs
		 	 	mslen = 0
		 	 	k = dmax

		 	 	for currentDeficit in range(dmax):
		 	 		while len(num_degs[k]) == 0:
		 	 			k -= 1
		 	 		target = num_degs[k].pop()
		 	 		G.add_edge(source, target)
		 	 		n -= 1
		 	 		if k > 1:
		 	 			modstubs[mslen] = (k-1,target)
		 	 			mslen += 1

		 	 	# Add back to the list any nonzero stubs that were removed
		 	 	for currentDeficit in range(mslen):
		 	 		(stubval, stubtarget) = modstubs[currentDeficit]
		 	 		num_degs[stubval].append(stubtarget)
		 	 		n += 1

		 G.name="havel_hakimi_graph %d nodes %d edges"%(G.order(),G.size())
		 return G
#endif



//		Aux::ShellList sl(seq);
//		index i = sl.size() - 1;
//
//		while (i > 0) {
//			while (! sl.isShellEmpty(i)) {
//				// get first vertex in shell
//				node v = sl.popVertexOfShell(i);
//				count deficit = i;
//				index shell = i;
//				std::set<node> tabu;
//
//				TRACE("node: ", v, ", shell: ", shell);
//
//				count nConnects = 0;
//				while (nConnects < deficit) {
//					// connect v to the other vertices
//					while (! sl.isShellEmpty(i)) {
//						node u = sl.popVertexOfShell(i);
//						if (tabu.count(u) > 0) {
//							// reinsert
//
//						}
//					}
//
//					// move connected vertices to appropriate shells
//
//					sl.forEachNodeInShell(shell, [&](node u) {
//						if (tabu.count(u) == 0 && nConnects < deficit) {
//							TRACE("try to add edge ", u, "/", v);
//							G.addEdge(u, v);
//							sl.decreaseShell(u);
//							++nConnects;
//							tabu.insert(u);
//						}
//					});
//					--shell;
//				}
//			}
//			--i;
//		}

//		-	Graph G2 = G;
//		-	while (G2.numberOfNodes() > 0) {
//		-		// TODO: main loop
//		+	for (int i = 0; i < sl.size(); i++) {
//		+		std::list<node>::iterator it = sl.getShelliterator(i);
//		+
//		+		sl.forEachNodeInShell(i, [&](node n) {
//		+			G.forNeighborsOf(n, [&](node m) {
//		+				if (sl.getCurrentShell(m) > i) {
//		+					sl.decreaseDegree(m);
//		+				}
//		+			});
//		+		});
//		 	}




		return G;
	}

}

} /* namespace NetworKit */

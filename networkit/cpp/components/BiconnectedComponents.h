/*
 * BiconnectedeComponents.h
 *
 * Created on: March 2018
 * 		 Author: Eugenio Angriman
 */

#ifndef BICONNECTEDCOMPONENTS_H_
#define BICONNECTEDCOMPONENTS_H_

#include "../base/Algorithm.h"
#include "../graph/Graph.h"

namespace NetworKit {

/*
 * @ingroup components
 * Determines the biconnected components of an undirected graph as defined in
 * Tarjan, Robert. Depth-First Search and Linear Graph Algorithms. SIAM J.
 * Comput. Vol 1, No. 2, June 1972.
 */
class BiconnectedComponents : public Algorithm {

public:
  /*
   * Creates BiconnectedComponents class for Graph @a G.
   *
   * @param G The graph.
   */
  BiconnectedComponents(const Graph &G);

  /*
   * This method determines the biconnected components for the graph given in
   * the constructor.
   */
  void run() override;

  /*
   *	Get the number of biconnected components.
   *
   *	@return The number of biconnected components.
   */
  count numberOfComponents();

  /*
   * Get the size of each component.
   *
   * @return Map from component index to size.
   */
  std::map<count, count> getComponentSizes();

  /*
   * Get the components vector.
   *
   * @return Vector of vectors, each component is stored as an (unordered) set
   * of nodes.
   */
  std::vector<std::vector<node>> getComponents();

protected:
  void init();
  void visit(node u);
  void newComponent(std::pair<node, node> e);

  const Graph &G;
  count n;
  count idx;
  count nComp;
  std::vector<count> level;
  std::vector<count> lowpt;
  std::vector<node> parent;
  std::vector<bool> visited;
  std::vector<bool> isRoot;
  std::vector<std::set<node>> componentsOfNode;
  std::map<count, count> componentSizes;
};

inline count BiconnectedComponents::numberOfComponents() {
  assureFinished();
  return nComp;
}

inline std::map<count, count> BiconnectedComponents::getComponentSizes() {
  assureFinished();
  return componentSizes;
}

inline std::vector<std::vector<node>> BiconnectedComponents::getComponents() {
  assureFinished();
  std::vector<std::vector<node>> result(nComp);

  node v = 0;
  for (auto components : componentsOfNode) {
    for (auto it = components.begin(); it != components.end(); ++it) {
      result[*it].push_back(v);
    }
    ++v;
  }
  return result;
}
} // namespace NetworKit
#endif // BICONNECTEDCOMPONENTS_H_

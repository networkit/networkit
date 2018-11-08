/*
* WeaklyConnectedComponents.h
*
*  Created on: June 20, 2017
*      Author: Eugenio Angriman
*/

#ifndef WEAKLYCONNECTEDCOMPONENTS_H_
#define WEAKLYCONNECTEDCOMPONENTS_H_

#include "../graph/Graph.h"
#include "../structures/Partition.h"
#include "../base/Algorithm.h"

namespace NetworKit {

    /**
    * @ingroup components
    * Determines the weakly connected components of a directed graph.
    */
    class WeaklyConnectedComponents : public Algorithm {
    public:
        /**
        * Create WeaklyConnectedComponents class for Graph @a G.
        *
        * @param G The graph.
        */
        WeaklyConnectedComponents(const Graph& G);

        /**
        * This method determines the weakly connected components for the graph
        * given in the constructor.
        */
        void run() override;

        /**
        * Get the number of weakly connected components.
        *
        * @return The number of weakly connected components.
        */
        count numberOfComponents();

        /**
        * Get the the component in which node @a u is.
        *
        * @param[in]	u	The index of the component where @a is located.
        */
        count componentOfNode(node u);

        /**
        * Get the map from component to size.
        *
        * @return A map that maps each component to its size.
        */
        std::map<index, count> getComponentSizes();

        /**
        * @return Vector of components, each stored as (unordered) set of nodes.
        */
        std::vector<std::vector<node>> getComponents();


    private:
        void updateComponent(
            index c, node w, std::queue<node>& q, bool inNeighbor
        );

        void init();

        // Pointer to the graph
        const Graph& G;

        // This vector associates each node to a component
        std::vector<index> components;

        // This map stores all components
        // <Key: component ID, Value: component size>
        std::map<index, count> compSize;
    };

    inline count WeaklyConnectedComponents::componentOfNode(node u) {
        assert (components[u] != none);
        assureFinished();
        return components[u];
    }

    inline count WeaklyConnectedComponents::numberOfComponents() {
        assureFinished();
        return compSize.size();
    }

    inline std::map<index, count> WeaklyConnectedComponents::getComponentSizes() {
       	assureFinished(); 
        return compSize;
    }

}


#endif /* WEAKLYCONNECTEDCOMPONENTS_H_ */

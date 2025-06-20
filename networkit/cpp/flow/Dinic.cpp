/*  Dinic.cpp
*
 *	Created on: 20.06.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <networkit/flow/Dinic.hpp>

namespace NetworKit{



void Dinic::run()
{
    buildResidual();
    edgeweight flow = 0;
    const edgeweight INF = std::numeric_limits<edgeweight>::max();

    // Main loop: while there is an augmenting path
    while (bfs()) {
        std::fill(ptr.begin(), ptr.end(), 0);
        // Send blocking flow
        while (edgeweight pushed = dfs(source, std::numeric_limits<edgeweight>::max()) {
            flow += pushed;
        }
    }
    maxFlowValue = flow;
    hasRun = true;
}
}

}
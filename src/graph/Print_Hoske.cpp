#include "Print_Hoske.h"

namespace NetworKit {

void print_hoske(std::ostream& os, const Graph& G) {
	G.forNodes([&] (node u) {
		os << u + 1 << "\n";
	});
	os << "#\n";
	G.forEdges([&] (node u, node v) {
		os << (u + 1) << " " << (v + 1) << "\n";
		os << (v + 1) << " " << (u + 1) << "\n";
	});
}

}

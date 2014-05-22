
%module NetworKit

%{
#include "../cpp/viz/Point.h"
#include "../cpp/graph/DirectedGraph.h"
#include "../cpp/graph/Graph.h"

#include "../cpp/graph/GraphGenerator.h"
#include "../cpp/generators/ChungLuGenerator.h"

#include "../cpp/properties/GraphProperties.h"
%}

%include exception.i

/* convert Python integer to C long */
%typemap(in) uint64_t {
	if (PyInt_Check($input)) {
		$1 = PyInt_AsLong($input);
	} else {
		SWIG_exception(SWIG_TypeError, "integer expected");
	}
}

/* use conversion rule in dispatch functions for overloaded methods */
%typemap(typecheck) uint64_t = int;

/* convert C long to Python integer */
%typemap(out) uint64_t {
    $result = PyInt_FromLong($1);
}

/* convert C++ std::pair to Python tuple */
%typemap(out) std::pair<uint64_t, uint64_t> {
	$result = PyTuple_New(2);
	PyTuple_SetItem($result, 0, PyInt_FromLong($1.first));
	PyTuple_SetItem($result, 1, PyInt_FromLong($1.second));
}


%include "../cpp/Globals.h"
%include "../cpp/viz/Point.h"
%include "../cpp/graph/IGraph.h"
%include "../cpp/graph/IDGraph.h"
%include "../cpp/graph/AbstractGraph.h"
%include "../cpp/graph/DirectedGraph.h"
%include "../cpp/graph/Graph.h"

%include "../cpp/structures/Partition.h"
%include "../cpp/graph/GraphGenerator.h"
%include "../cpp/generators/ChungLuGenerator.h"

%include "../cpp/properties/GraphProperties.h"

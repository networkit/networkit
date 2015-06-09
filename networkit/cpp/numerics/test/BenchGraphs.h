/*
 * BenchGraphs.h
 *
 *  Created on: May 27, 2015
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_CPP_NUMERICS_TEST_BENCHGRAPHS_H_
#define NETWORKIT_CPP_NUMERICS_TEST_BENCHGRAPHS_H_

#include <vector>
#include <string>

using namespace std;

const vector<string> GRIDS {
	"instances/grid/Laplace_4x4.graph",
	"instances/grid/Laplace_8x8.graph",
	"instances/grid/Laplace_16x16.graph",
	"instances/grid/Laplace_32x32.graph",
	"instances/grid/Laplace_64x64.graph",
	"instances/grid/Laplace_128x128.graph",
	"instances/grid/Laplace_256x256.graph",
	"instances/grid/Laplace_512x512.graph",
	"instances/grid/Laplace_1024x1024.graph"
};

const vector<string> BARABASI {
	"instances/barabasi/25_att_4_unweighted.graph",
	"instances/barabasi/100_att_4_unweighted.graph",
	"instances/barabasi/500_att_4_unweighted.graph",
	"instances/barabasi/1000_att_4_unweighted.graph",
	"instances/barabasi/5000_att_4_unweighted.graph",
	"instances/barabasi/10000_att_4_unweighted.graph",
	"instances/barabasi/50000_att_4_unweighted.graph",
};

const vector<string> LAPLACE { // https://www.paralution.com/downloads/Laplace/
	"instances/laplace/Laplace_5x5.mtx",
	"instances/laplace/Laplace_20x20.mtx",
	"instances/laplace/Laplace_40x40.mtx",
	"instances/laplace/Laplace_80x80.mtx",
	"instances/laplace/Laplace_200x200.mtx",
	"instances/laplace/Laplace_300x300.mtx",
	"instances/laplace/Laplace_400x400.mtx",
	"instances/laplace/Laplace_500x500.mtx",
	"instances/laplace/Laplace_600x600.mtx",
	"instances/laplace/Laplace_700x700.mtx",
	"instances/laplace/Laplace_800x800.mtx",
	"instances/laplace/Laplace_900x900.mtx",
	"instances/laplace/Laplace_1000x1000.mtx"
};

const vector<string> PITZ_DAILY { // https://www.paralution.com/downloads/benchmarks/pitzDaily_16M_MPI1.tar.gz
	"instances/pitzDaily_16M_MPI1/pitzDaily_mat_1.csr",
	"instances/pitzDaily_16M_MPI1/pitzDaily_mat_400.csr"
};

const vector<string> DIMACS_NUMERICS { // http://www.cc.gatech.edu/dimacs10/archive/numerical.shtml
	"instances/dimacs_numerics/333Sp.graph",
	"instances/dimacs_numerics/AS365.graph",
	"instances/dimacs_numerics/M6.graph",
	"instances/dimacs_numerics/NACA0015.graph",
	"instances/dimacs_numerics/NLR.graph"
};


const vector<string> DIMACS_SPARSE { // http://www.cc.gatech.edu/dimacs10/archive/matrix.shtml
	"instances/dimacs_sparse_matrices/af_shell9.graph",
	"instances/dimacs_sparse_matrices/af_shell10.graph",
	"instances/dimacs_sparse_matrices/audikw1.graph",
	"instances/dimacs_sparse_matrices/cage15.graph",
	"instances/dimacs_sparse_matrices/ecology1.graph",
	"instances/dimacs_sparse_matrices/ecology2.mtx",
	"instances/dimacs_sparse_matrices/G3_circuit.graph",
	"instances/dimacs_sparse_matrices/kkt_power.graph",
	"instances/dimacs_sparse_matrices/ldoor.graph",
	"instances/dimacs_sparse_matrices/nlpkkt120.graph",
	"instances/dimacs_sparse_matrices/nlpkkt160.graph", // 1,75 GB
	"instances/dimacs_sparse_matrices/nlpkkt200.graph", // 3,61 GB
	"instances/dimacs_sparse_matrices/thermal2.graph"
};

const vector<string> WALSHAW { // http://staffweb.cms.gre.ac.uk/~wc06/partition/
	"instances/walshaw/3elt.graph",
	"instances/walshaw/4elt.graph",
	"instances/walshaw/144.graph",
	"instances/walshaw/598a.graph",
	"instances/walshaw/add20.graph",
	"instances/walshaw/add32.graph",
	"instances/walshaw/auto.graph",
	"instances/walshaw/bcsstk29.graph",
	"instances/walshaw/bcsstk30.graph",
	"instances/walshaw/bcsstk31.graph",
	"instances/walshaw/bcsstk32.graph",
	"instances/walshaw/bcsstk33.graph",
	"instances/walshaw/brack2.graph",
	"instances/walshaw/crack.graph",
	"instances/walshaw/cs4.graph",
	"instances/walshaw/cti.graph",
	"instances/walshaw/data.graph",
	"instances/walshaw/fe_4elt2.graph",
	"instances/walshaw/fe_body.graph",
	"instances/walshaw/fe_ocean.graph",
	"instances/walshaw/fe_pwt.graph",
	"instances/walshaw/fe_rotor.graph",
	"instances/walshaw/fe_sphere.graph",
	"instances/walshaw/fe_tooth.graph",
	"instances/walshaw/finan512.graph",
	"instances/walshaw/m14b.graph",
	"instances/walshaw/memplus.graph",
	"instances/walshaw/t60k.graph",
	"instances/walshaw/uk.graph",
	"instances/walshaw/vibrobox.graph",
	"instances/walshaw/wave.graph",
	"instances/walshaw/whitaker3.graph",
	"instances/walshaw/wing_nodal.graph",
	"instances/walshaw/wing.graph"
};



#endif /* NETWORKIT_CPP_NUMERICS_TEST_BENCHGRAPHS_H_ */

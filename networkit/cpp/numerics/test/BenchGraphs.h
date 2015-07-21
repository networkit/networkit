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

enum Type {MATRIX, NETWORK, LAPLACIAN, SDD, UNKNOWN};

struct Instance {
	string path;
	Type type;
};

const vector<Instance> GRIDS {
	{"instances/grid/Laplace_4x4.graph", NETWORK},
	{"instances/grid/Laplace_8x8.graph", NETWORK},
	{"instances/grid/Laplace_16x16.graph", NETWORK},
	{"instances/grid/Laplace_32x32.graph", NETWORK},
	{"instances/grid/Laplace_64x64.graph", NETWORK},
	{"instances/grid/Laplace_128x128.graph", NETWORK},
	{"instances/grid/Laplace_256x256.graph", NETWORK},
	{"instances/grid/Laplace_512x512.graph", NETWORK},
	{"instances/grid/Laplace_1024x1024.graph", NETWORK}
};

const vector<Instance> BARABASI {
	{"instances/barabasi/25_att_4_unweighted.graph", NETWORK},
	{"instances/barabasi/100_att_4_unweighted.graph", NETWORK},
	{"instances/barabasi/500_att_4_unweighted.graph", NETWORK},
	{"instances/barabasi/1000_att_4_unweighted.graph", NETWORK},
	{"instances/barabasi/5000_att_4_unweighted.graph", NETWORK},
	{"instances/barabasi/10000_att_4_unweighted.graph", NETWORK},
	{"instances/barabasi/50000_att_4_unweighted.graph", NETWORK}
};

const vector<Instance> LAPLACE { // https://www.paralution.com/downloads/Laplace/
	{"instances/laplace/Laplace_5x5.mtx", MATRIX},
	{"instances/laplace/Laplace_20x20.mtx", MATRIX},
	{"instances/laplace/Laplace_40x40.mtx", MATRIX},
	{"instances/laplace/Laplace_80x80.mtx", MATRIX},
	{"instances/laplace/Laplace_200x200.mtx", MATRIX},
	{"instances/laplace/Laplace_300x300.mtx", MATRIX},
	{"instances/laplace/Laplace_400x400.mtx", MATRIX},
	{"instances/laplace/Laplace_500x500.mtx", MATRIX},
	{"instances/laplace/Laplace_600x600.mtx", MATRIX},
	{"instances/laplace/Laplace_700x700.mtx", MATRIX},
	{"instances/laplace/Laplace_800x800.mtx", MATRIX},
	{"instances/laplace/Laplace_900x900.mtx", MATRIX},
	{"instances/laplace/Laplace_1000x1000.mtx", MATRIX}
};

const vector<Instance> PITZ_DAILY { // https://www.paralution.com/downloads/benchmarks/pitzDaily_16M_MPI1.tar.gz
	{"instances/pitzDaily_16M_MPI1/pitzDaily_mat_1.csr", MATRIX},
	{"instances/pitzDaily_16M_MPI1/pitzDaily_mat_400.csr", MATRIX}
};

const vector<Instance> DIMACS_NUMERICS { // http://www.cc.gatech.edu/dimacs10/archive/numerical.shtml
	{"instances/dimacs_numerics/333Sp.graph", MATRIX},
	{"instances/dimacs_numerics/AS365.graph", MATRIX},
	{"instances/dimacs_numerics/M6.graph", MATRIX},
	{"instances/dimacs_numerics/NACA0015.graph", MATRIX},
	{"instances/dimacs_numerics/NLR.graph", MATRIX}
};


const vector<Instance> DIMACS_SPARSE { // http://www.cc.gatech.edu/dimacs10/archive/matrix.shtml
//	{"instances/dimacs_sparse_matrices/af_shell9.mtx", MATRIX}, // no sdd matrix
//	{"instances/dimacs_sparse_matrices/af_shell10.mtx", MATRIX}, // no sdd matrix
//	{"instances/dimacs_sparse_matrices/audikw_1.mtx", MATRIX}, // no sdd matrix
//	{"instances/dimacs_sparse_matrices/cage15.mtx", MATRIX}, // no sdd matrix
	{"instances/dimacs_sparse_matrices/ecology1.mtx", MATRIX},
	{"instances/dimacs_sparse_matrices/ecology2.mtx", MATRIX},
//	{"instances/dimacs_sparse_matrices/G3_circuit.mtx", MATRIX},
//	{"instances/dimacs_sparse_matrices/kkt_power.mtx", MATRIX}, // no sdd matrix
//	{"instances/dimacs_sparse_matrices/ldoor.graph", MATRIX},
//	{"instances/dimacs_sparse_matrices/nlpkkt120.mtx", MATRIX}, // no sdd matrix
//	{"instances/dimacs_sparse_matrices/nlpkkt160.mtx", MATRIX}, // 1,75 GB, no sdd matrix
//	{"instances/dimacs_sparse_matrices/nlpkkt200.graph", MATRIX}, // 3,61 GB
//	{"instances/dimacs_sparse_matrices/thermal2.mtx", MATRIX} // no sdd matrix
};

const vector<Instance> DIMACS_CLUSTERING { // http://www.cc.gatech.edu/dimacs10/archive/clustering.shtml
	{"instances/dimacs_clustering/as-22july06.graph", NETWORK},
	{"instances/dimacs_clustering/cnr-2000.graph", NETWORK},
	{"instances/dimacs_clustering/in-2004.graph", NETWORK},
	{"instances/dimacs_clustering/PGPgiantcompo.graph", NETWORK},
	{"instances/dimacs_clustering/preferentialAttachment.graph", NETWORK},
	{"instances/dimacs_clustering/smallworld.graph", NETWORK},
	{"instances/dimacs_clustering/uk-2002.graph", NETWORK}
};

const vector<Instance> DIMACS_STREET_NETWORKS { // http://www.cc.gatech.edu/dimacs10/archive/streets.shtml
	{"instances/dimacs_street_networks/belgium.osm.graph", NETWORK},
	{"instances/dimacs_street_networks/europe.osm.graph", NETWORK},
	{"instances/dimacs_street_networks/great-britain.osm.graph", NETWORK},
	{"instances/dimacs_street_networks/luxembourg.osm.graph", NETWORK}
};

const vector<Instance> FACEBOOK100 {
//	{"instances/facebook100/American75.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Amherst41.mat.txt.metis", NETWORK},
	{"instances/facebook100/Auburn71.mat.txt.metis", NETWORK},
//	{"instances/facebook100/BC17.mat.txt.metis", NETWORK},
//	{"instances/facebook100/BU10.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Baylor93.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Berkeley13.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Bingham82.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Bowdoin47.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Brandeis99.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Brown11.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Bucknell39.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Cal65.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Caltech36.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Carnegie49.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Colgate88.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Columbia2.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Cornell5.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Dartmouth6.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Duke14.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Emory27.mat.txt.metis", NETWORK},
//	{"instances/facebook100/FSU53.mat.txt.metis", NETWORK},
//	{"instances/facebook100/GWU54.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Georgetown15.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Hamilton46.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Harvard1.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Haverford76.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Howard90.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Indiana69.mat.txt.metis", NETWORK},
//	{"instances/facebook100/JMU79.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Lehigh96.mat.txt.metis", NETWORK},
//	{"instances/facebook100/MIT8.mat.txt.metis", NETWORK},
//	{"instances/facebook100/MSU24.mat.txt.metis", NETWORK},
//	{"instances/facebook100/MU78.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Maine59.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Maryland58.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Mich67.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Michigan23.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Middlebury45.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Mississippi66.mat.txt.metis", NETWORK},
//	{"instances/facebook100/NYU9.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Northeastern19.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Northwestern25.mat.txt.metis", NETWORK},
//	{"instances/facebook100/NotreDame57.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Oberlin44.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Oklahoma97.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Penn94.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Pepperdine86.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Princeton12.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Reed98.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Rice31.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Rochester38.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Rutgers89.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Santa74.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Simmons81.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Smith60.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Stanford3.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Swarthmore42.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Syracuse56.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Temple83.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Tennessee95.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Texas80.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Texas84.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Trinity100.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Tufts18.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Tulane29.mat.txt.metis", NETWORK},
//	{"instances/facebook100/UC33.mat.txt.metis", NETWORK},
//	{"instances/facebook100/UC61.mat.txt.metis", NETWORK},
//	{"instances/facebook100/UC64.mat.txt.metis", NETWORK},
//	{"instances/facebook100/UCF52.mat.txt.metis", NETWORK},
//	{"instances/facebook100/UCLA26.mat.txt.metis", NETWORK},
//	{"instances/facebook100/UCSB37.mat.txt.metis", NETWORK},
//	{"instances/facebook100/UCSC68.mat.txt.metis", NETWORK},
//	{"instances/facebook100/UCSD34.mat.txt.metis", NETWORK},
//	{"instances/facebook100/UChicago30.mat.txt.metis", NETWORK},
//	{"instances/facebook100/UConn91.mat.txt.metis", NETWORK},
//	{"instances/facebook100/UF21.mat.txt.metis", NETWORK},
//	{"instances/facebook100/UGA50.mat.txt.metis", NETWORK},
//	{"instances/facebook100/UIllinois20.mat.txt.metis", NETWORK},
//	{"instances/facebook100/UMass92.mat.txt.metis", NETWORK},
//	{"instances/facebook100/UNC28.mat.txt.metis", NETWORK},
//	{"instances/facebook100/UPenn7.mat.txt.metis", NETWORK},
//	{"instances/facebook100/USC35.mat.txt.metis", NETWORK},
//	{"instances/facebook100/USF51.mat.txt.metis", NETWORK},
//	{"instances/facebook100/USFCA72.mat.txt.metis", NETWORK},
//	{"instances/facebook100/UVA16.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Vanderbilt48.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Vassar85.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Vermont70.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Villanova62.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Virginia63.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Wake73.mat.txt.metis", NETWORK},
//	{"instances/facebook100/WashU32.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Wellesley22.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Wesleyan43.mat.txt.metis", NETWORK},
//	{"instances/facebook100/William77.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Williams40.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Wisconsin87.mat.txt.metis", NETWORK},
//	{"instances/facebook100/Yale4.mat.txt.metis", NETWORK},
};


const vector<Instance> WALSHAW { // http://staffweb.cms.gre.ac.uk/~wc06/partition/
	{"instances/walshaw/3elt.graph", NETWORK},
	{"instances/walshaw/4elt.graph", NETWORK},
	{"instances/walshaw/144.graph", NETWORK},
	{"instances/walshaw/598a.graph", NETWORK},
	{"instances/walshaw/add20.graph", NETWORK},
	{"instances/walshaw/add32.graph", NETWORK},
	{"instances/walshaw/auto.graph", NETWORK},
	{"instances/walshaw/bcsstk29.graph", NETWORK},
	{"instances/walshaw/bcsstk30.graph", NETWORK},
	{"instances/walshaw/bcsstk31.graph", NETWORK},
	{"instances/walshaw/bcsstk32.graph", NETWORK},
	{"instances/walshaw/bcsstk33.graph", NETWORK},
	{"instances/walshaw/brack2.graph", NETWORK},
	{"instances/walshaw/crack.graph", NETWORK},
	{"instances/walshaw/cs4.graph", NETWORK},
	{"instances/walshaw/cti.graph", NETWORK},
	{"instances/walshaw/data.graph", NETWORK},
	{"instances/walshaw/fe_4elt2.graph", NETWORK},
	{"instances/walshaw/fe_body.graph", NETWORK},
	{"instances/walshaw/fe_ocean.graph", NETWORK},
	{"instances/walshaw/fe_pwt.graph", NETWORK},
	{"instances/walshaw/fe_rotor.graph", NETWORK},
	{"instances/walshaw/fe_sphere.graph", NETWORK},
	{"instances/walshaw/fe_tooth.graph", NETWORK},
	{"instances/walshaw/finan512.graph", NETWORK},
	{"instances/walshaw/m14b.graph", NETWORK},
	{"instances/walshaw/memplus.graph", NETWORK},
	{"instances/walshaw/t60k.graph", NETWORK},
	{"instances/walshaw/uk.graph", NETWORK},
	{"instances/walshaw/vibrobox.graph", NETWORK},
	{"instances/walshaw/wave.graph", NETWORK},
	{"instances/walshaw/whitaker3.graph", NETWORK},
	{"instances/walshaw/wing_nodal.graph", NETWORK},
	{"instances/walshaw/wing.graph", NETWORK}
};






#endif /* NETWORKIT_CPP_NUMERICS_TEST_BENCHGRAPHS_H */

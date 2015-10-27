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
	bool isConnected;
};

const vector<Instance> GRIDS {
	{"instances/grid/Laplace_4x4.graph", NETWORK, true},
	{"instances/grid/Laplace_8x8.graph", NETWORK, true},
	{"instances/grid/Laplace_16x16.graph", NETWORK, true},
	{"instances/grid/Laplace_32x32.graph", NETWORK, true},
	{"instances/grid/Laplace_64x64.graph", NETWORK, true},
	{"instances/grid/Laplace_128x128.graph", NETWORK, true},
	{"instances/grid/Laplace_256x256.graph", NETWORK, true},
	{"instances/grid/Laplace_512x512.graph", NETWORK, true},
	{"instances/grid/Laplace_1024x1024.graph", NETWORK, true}
};

const vector<Instance> BARABASI {
	{"instances/barabasi/25_att_4_unweighted.graph", NETWORK, true},
	{"instances/barabasi/100_att_4_unweighted.graph", NETWORK, true},
	{"instances/barabasi/500_att_4_unweighted.graph", NETWORK, true},
	{"instances/barabasi/1000_att_4_unweighted.graph", NETWORK, true},
	{"instances/barabasi/5000_att_4_unweighted.graph", NETWORK, true},
	{"instances/barabasi/10000_att_4_unweighted.graph", NETWORK, true},
	{"instances/barabasi/50000_att_4_unweighted.graph", NETWORK, true}
};

const vector<Instance> LAPLACE { // https://www.paralution.com/downloads/Laplace/
	{"instances/laplace/Laplace_5x5.mtx", MATRIX, true},
	{"instances/laplace/Laplace_20x20.mtx", MATRIX, true},
	{"instances/laplace/Laplace_40x40.mtx", MATRIX, true},
	{"instances/laplace/Laplace_80x80.mtx", MATRIX, true},
	{"instances/laplace/Laplace_200x200.mtx", MATRIX, true},
	{"instances/laplace/Laplace_300x300.mtx", MATRIX, true},
	{"instances/laplace/Laplace_400x400.mtx", MATRIX, true},
	{"instances/laplace/Laplace_500x500.mtx", MATRIX, true},
	{"instances/laplace/Laplace_600x600.mtx", MATRIX, true},
	{"instances/laplace/Laplace_700x700.mtx", MATRIX, true},
	{"instances/laplace/Laplace_800x800.mtx", MATRIX, true},
	{"instances/laplace/Laplace_900x900.mtx", MATRIX, true},
	{"instances/laplace/Laplace_1000x1000.mtx", MATRIX, true}
};

const vector<Instance> PITZ_DAILY { // https://www.paralution.com/downloads/benchmarks/pitzDaily_16M_MPI1.tar.gz
	{"instances/pitzDaily_16M_MPI1/pitzDaily_mat_1.csr", MATRIX, true},
	{"instances/pitzDaily_16M_MPI1/pitzDaily_mat_400.csr", MATRIX, true}
};

const vector<Instance> DIMACS_NUMERICS { // http://www.cc.gatech.edu/dimacs10/archive/numerical.shtml
	{"instances/dimacs_numerics/333Sp.graph", MATRIX, true},
	{"instances/dimacs_numerics/AS365.graph", MATRIX, true},
	{"instances/dimacs_numerics/M6.graph", MATRIX, true},
	{"instances/dimacs_numerics/NACA0015.graph", MATRIX, true},
	{"instances/dimacs_numerics/NLR.graph", MATRIX, true}
};


const vector<Instance> DIMACS_SPARSE { // http://www.cc.gatech.edu/dimacs10/archive/MATRIX, true.shtml
//	{"instances/dimacs_sparse_matrices/af_shell9.mtx", MATRIX, true}, // no sdd MATRIX, true
//	{"instances/dimacs_sparse_matrices/af_shell10.mtx", MATRIX, true}, // no sdd MATRIX, true
//	{"instances/dimacs_sparse_matrices/audikw_1.mtx", MATRIX, true}, // no sdd MATRIX, true
//	{"instances/dimacs_sparse_matrices/cage15.mtx", MATRIX, true}, // no sdd MATRIX, true
	{"instances/dimacs_sparse_matrices/ecology1.mtx", MATRIX, true},
	{"instances/dimacs_sparse_matrices/ecology2.mtx", MATRIX, true},
//	{"instances/dimacs_sparse_matrices/G3_circuit.mtx", MATRIX, true},
//	{"instances/dimacs_sparse_matrices/kkt_power.mtx", MATRIX, true}, // no sdd MATRIX, true
//	{"instances/dimacs_sparse_matrices/ldoor.graph", MATRIX, true},
//	{"instances/dimacs_sparse_matrices/nlpkkt120.mtx", MATRIX, true}, // no sdd MATRIX, true
//	{"instances/dimacs_sparse_matrices/nlpkkt160.mtx", MATRIX, true}, // 1,75 GB, no sdd MATRIX, true
//	{"instances/dimacs_sparse_matrices/nlpkkt200.graph", MATRIX, true}, // 3,61 GB
//	{"instances/dimacs_sparse_matrices/thermal2.mtx", MATRIX, true} // no sdd MATRIX, true
};

const vector<Instance> DIMACS_CITATION {
	{"instances/dimacs_citation_networks/citationCiteseer.graph", NETWORK, true},
	{"instances/dimacs_citation_networks/coAuthorsCiteseer.graph", NETWORK, true},
	{"instances/dimacs_citation_networks/coAuthorsDBLP.graph", NETWORK, true},
	{"instances/dimacs_citation_networks/coPapersCiteseer.graph", NETWORK, true},
	{"instances/dimacs_citation_networks/coPapersDBLP.graph", NETWORK, true},	
};

const vector<Instance> DIMACS_CLUSTERING { // http://www.cc.gatech.edu/dimacs10/archive/clustering.shtml
	{"instances/dimacs_clustering/as-22july06.graph", NETWORK, true},
	{"instances/dimacs_clustering/astro-ph.graph", NETWORK, true},
	{"instances/dimacs_clustering/cnr-2000.graph", NETWORK, true},
	{"instances/dimacs_clustering/hep-th.graph", NETWORK, true},
	{"instances/dimacs_clustering/in-2004.graph", NETWORK, true},
	{"instances/dimacs_clustering/PGPgiantcompo.graph", NETWORK, true},
	{"instances/dimacs_clustering/preferentialAttachment.graph", NETWORK, true},
	{"instances/dimacs_clustering/smallworld.graph", NETWORK, true}
//	{"instances/dimacs_clustering/uk-2002.graph", NETWORK, true}
};

const vector<Instance> DIMACS_STREET_NETWORKS { // http://www.cc.gatech.edu/dimacs10/archive/streets.shtml
	{"instances/dimacs_street_networks/asia.osm.graph", NETWORK, true},
	{"instances/dimacs_street_networks/belgium.osm.graph", NETWORK, true},
	{"instances/dimacs_street_networks/europe.osm.graph", NETWORK, true},
	{"instances/dimacs_street_networks/germany.osm.graph", NETWORK, true},
	{"instances/dimacs_street_networks/great-britain.osm.graph", NETWORK, true},
	{"instances/dimacs_street_networks/luxembourg.osm.graph", NETWORK, true},
	{"instances/dimacs_street_networks/netherlands.osm.graph", NETWORK, true}
};

const vector<Instance> FACEBOOK100 {
	{"instances/facebook100/American75.mat.txt.metis", NETWORK, false},
	{"instances/facebook100/Amherst41.mat.txt.metis", NETWORK, true},
	{"instances/facebook100/Auburn71.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/BC17.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/BU10.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Baylor93.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Berkeley13.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Bingham82.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Bowdoin47.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Brandeis99.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Brown11.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Bucknell39.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Cal65.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Caltech36.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Carnegie49.mat.txt.metis", NETWORK, true},
	{"instances/facebook100/Colgate88.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Columbia2.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Cornell5.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Dartmouth6.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Duke14.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Emory27.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/FSU53.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/GWU54.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Georgetown15.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Hamilton46.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Harvard1.mat.txt.metis", NETWORK, true},
	{"instances/facebook100/Haverford76.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Howard90.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Indiana69.mat.txt.metis", NETWORK, true},
	{"instances/facebook100/JMU79.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Lehigh96.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/MIT8.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/MSU24.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/MU78.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Maine59.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Maryland58.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Mich67.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Michigan23.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Middlebury45.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Mississippi66.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/NYU9.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Northeastern19.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Northwestern25.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/NotreDame57.mat.txt.metis", NETWORK, true},
	{"instances/facebook100/Oberlin44.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Oklahoma97.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Penn94.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Pepperdine86.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Princeton12.mat.txt.metis", NETWORK, true},
	{"instances/facebook100/Reed98.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Rice31.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Rochester38.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Rutgers89.mat.txt.metis", NETWORK, true},
	{"instances/facebook100/Santa74.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Simmons81.mat.txt.metis", NETWORK, true},
	{"instances/facebook100/Smith60.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Stanford3.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Swarthmore42.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Syracuse56.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Temple83.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Tennessee95.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Texas80.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Texas84.mat.txt.metis", NETWORK, true},
	{"instances/facebook100/Trinity100.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Tufts18.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Tulane29.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/UC33.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/UC61.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/UC64.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/UCF52.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/UCLA26.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/UCSB37.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/UCSC68.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/UCSD34.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/UChicago30.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/UConn91.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/UF21.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/UGA50.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/UIllinois20.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/UMass92.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/UNC28.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/UPenn7.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/USC35.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/USF51.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/USFCA72.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/UVA16.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Vanderbilt48.mat.txt.metis", NETWORK, true},
	{"instances/facebook100/Vassar85.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Vermont70.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Villanova62.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Virginia63.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Wake73.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/WashU32.mat.txt.metis", NETWORK, true},
	{"instances/facebook100/Wellesley22.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Wesleyan43.mat.txt.metis", NETWORK, true},
	{"instances/facebook100/William77.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Williams40.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Wisconsin87.mat.txt.metis", NETWORK, true},
//	{"instances/facebook100/Yale4.mat.txt.metis", NETWORK, true},
};


const vector<Instance> WALSHAW { // http://staffweb.cms.gre.ac.uk/~wc06/partition/
	{"instances/walshaw/3elt.graph", NETWORK, true},
	{"instances/walshaw/4elt.graph", NETWORK, true},
	{"instances/walshaw/144.graph", NETWORK, true},
	{"instances/walshaw/598a.graph", NETWORK, true},
	{"instances/walshaw/add20.graph", NETWORK, true},
	{"instances/walshaw/add32.graph", NETWORK, true},
	{"instances/walshaw/auto.graph", NETWORK, true},
	{"instances/walshaw/bcsstk29.graph", NETWORK, false},
	{"instances/walshaw/bcsstk30.graph", NETWORK, false},
	{"instances/walshaw/bcsstk31.graph", NETWORK, false},
	{"instances/walshaw/bcsstk32.graph", NETWORK, true},
	{"instances/walshaw/bcsstk33.graph", NETWORK, true},
	{"instances/walshaw/brack2.graph", NETWORK, true},
	{"instances/walshaw/crack.graph", NETWORK, true},
	{"instances/walshaw/cs4.graph", NETWORK, true},
	{"instances/walshaw/cti.graph", NETWORK, true},
	{"instances/walshaw/data.graph", NETWORK, true},
	{"instances/walshaw/fe_4elt2.graph", NETWORK, true},
	{"instances/walshaw/fe_body.graph", NETWORK, false},
	{"instances/walshaw/fe_ocean.graph", NETWORK, true},
	{"instances/walshaw/fe_pwt.graph", NETWORK, false},
	{"instances/walshaw/fe_rotor.graph", NETWORK, true},
	{"instances/walshaw/fe_sphere.graph", NETWORK, true},
	{"instances/walshaw/fe_tooth.graph", NETWORK, true},
	{"instances/walshaw/finan512.graph", NETWORK, true},
	{"instances/walshaw/m14b.graph", NETWORK, true},
	{"instances/walshaw/memplus.graph", NETWORK, true},
	{"instances/walshaw/t60k.graph", NETWORK, true},
	{"instances/walshaw/uk.graph", NETWORK, true},
	{"instances/walshaw/vibrobox.graph", NETWORK, true},
	{"instances/walshaw/wave.graph", NETWORK, true},
	{"instances/walshaw/whitaker3.graph", NETWORK, true},
	{"instances/walshaw/wing_nodal.graph", NETWORK, true},
	{"instances/walshaw/wing.graph", NETWORK, true}
};

const vector<Instance> SNAP {
//	{"Brightkite_edges.snap", NETWORK, true},
//	{"CA-AstroPh.snap", NETWORK, true},
//	{"CA-CondMat.snap", NETWORK, true},
	{"/Users/Michael/Downloads/SNAP/CA-GrQc.graph", NETWORK, false},
//	{"/Users/Michael/Downloads/SNAP/CA-HepPh.graph", NETWORK, true},
//	{"/Users/Michael/Downloads/SNAP/CA-HepTh.sna.graph", NETWORK, true},
//	{"Email-Enron.snap", NETWORK, true},
//	{"Gowalla_edges.snap", NETWORK, true},
//	{"as-skitter.snap", NETWORK, true},
	{"/Users/Michael/Downloads/SNAP/com-amazon.ungraph.graph", NETWORK, false},
//	{"com-orkut.ungraph.snap", NETWORK, true},
//	{"com-youtube.ungraph.snap", NETWORK, true},
//	{"roadNet-CA.snap", NETWORK, true},
//	{"roadNet-PA.snap", NETWORK, true},
//	{"roadNet-TX.snap", NETWORK, true},
};





#endif /* NETWORKIT_CPP_NUMERICS_TEST_BENCHGRAPHS_H */

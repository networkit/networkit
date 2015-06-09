/**
 * Set current benchmarks here
 */

const vector<string> GRIDS = {
  "grid/Laplace_5x5.mtx",
  "grid/Laplace_10x10.mtx",
  "grid/Laplace_20x20.mtx",
  "grid/Laplace_30x30.mtx",
  "grid/Laplace_40x40.mtx",
  "grid/Laplace_50x50.mtx",
  "grid/Laplace_60x60.mtx",
  "grid/Laplace_70x70.mtx",
  "grid/Laplace_80x80.mtx",
  "grid/Laplace_90x90.mtx",
  "grid/Laplace_100x100.mtx",
  "grid/Laplace_200x200.mtx",
  "grid/Laplace_250x250.graph",
  "grid/Laplace_300x300.mtx",
  "grid/Laplace_350x350.graph",
  "grid/Laplace_400x400.mtx",
  "grid/Laplace_450x450.graph",
  "grid/Laplace_500x500.mtx",
  "grid/Laplace_550x550.graph",
  "grid/Laplace_600x600.mtx",
  "grid/Laplace_650x650.graph",
  "grid/Laplace_700x700.mtx",
  "grid/Laplace_750x750.graph",
  "grid/Laplace_800x800.mtx",
  "grid/Laplace_850x850.graph",
  "grid/Laplace_900x900.mtx",
  "grid/Laplace_950x950.graph",
  "grid/Laplace_1000x1000.mtx",
};
const vector<string> GRIDS_3D = {
  "grid3/Laplace_10x10x10.graph",
  "grid3/Laplace_20x20x20.graph",
  "grid3/Laplace_30x30x30.graph",
  "grid3/Laplace_40x40x40.graph",
  "grid3/Laplace_50x50x50.graph",
  "grid3/Laplace_60x60x60.graph",
  "grid3/Laplace_70x70x70.graph",
  "grid3/Laplace_75x75x75.graph",
  "grid3/Laplace_80x80x80.graph",
  "grid3/Laplace_85x85x85.graph",
  "grid3/Laplace_90x90x90.graph",
  "grid3/Laplace_95x95x95.graph",
  "grid3/Laplace_100x100x100.graph",
};
const vector<string> BARABASI = {
  "barabasi/100_att_4.graph",
  "barabasi/500_att_4.graph",
  "barabasi/1000_att_4.graph",
  "barabasi/5000_att_4.graph",
  "barabasi/10000_att_4.graph",
  "barabasi/25000_att_4.graph",
  "barabasi/50000_att_4.graph",
  "barabasi/100000_att_4.graph",
  "barabasi/150000_att_4.graph",
  "barabasi/200000_att_4.graph",
  "barabasi/250000_att_4.graph",
  "barabasi/300000_att_4.graph",
  "barabasi/350000_att_4.graph",
  "barabasi/400000_att_4.graph",
  "barabasi/450000_att_4.graph",
  "barabasi/500000_att_4.graph",
  "barabasi/600000_att_4.graph",
  "barabasi/700000_att_4.graph",
  "barabasi/800000_att_4.graph",
  "barabasi/900000_att_4.graph",
  "barabasi/1000000_att_4.graph",
};


const vector<string> BENCH_GRAPHS = {
  "lesmis.graph",
  "power.graph",
//  //    {"luxembourg.osm.graph", "http://www.cc.gatech.edu/dimacs10/archive/streets.shtml"},
//  //  {"shallow-water1.mtx", "http://www.cise.ufl.edu/research/sparse/matrices/MaxPlanck/shallow_water1.html"},
    "caidaRouterLevel.graph",
//  //  {"cnr-2000.graph", "http://www.cc.gatech.edu/dimacs10/archive/clustering.shtml"},
//  //  {"parabolic-fem.mtx", "http://www.cise.ufl.edu/research/sparse/matrices/Wissgott/parabolic_fem.html"},
//  //  {"NACA0015.graph", "http://www.cc.gatech.edu/dimacs10/archive/numerical.shtml"},
//  //  {"amazon-ratings.tsv", "http://konect.uni-koblenz.de/networks/amazon-ratings"},
//  //  {"hugetrace-00000.graph", "http://www.cc.gatech.edu/dimacs10/archive/dyn-frames.shtml"},
};

const vector<Benchmark> BENCHS = {
//  {
//    "DEBUG",
//    BenchType::DEBUG,
//    5,
//    1e-4,
//    123456,
//    BENCH_GRAPHS,
//    {}
//  },
//  {
//    "tree",
//    BenchType::TREE_COMPLEXITY,
//    1,
//    1e-2,
//    123456,
//    {"grid/Laplace_100x100.mtx", "grid/Laplace_100x100_weighted.graph",
//     //"luxembourg.osm.graph",
//     //"barabasi/25000_att_4_unweighted.graph", "barabasi/25000_att_4.graph",
//     //"caidaRouterLevel.graph"
//    },
//    {
//    }
//  },
//  {
//    "tree_flow",
//    BenchType::TREE_COMPLEXITY,
//    1,
//    1e-2,
//    123456,
//    {"grid/Laplace_1000x1000.mtx", //"grid/Laplace_500x500.mtx",
//     //"barabasi/25000_att_4.graph", "barabasi/200000_att_4.graph",
//    },
//    {
//    }
//  },
//  {
//    "grid3_tree",
//    BenchType::TREE_COMPLEXITY,
//    1,
//    1e-2,
//    123456,
//    {"grid3/Laplace_50x50x50.graph"},
//    {
//    }
//  },
//  {
//    "barabasi_tree",
//    BenchType::TREE_COMPLEXITY,
//    1,
//    1e-2,
//    123456,
//    {"barabasi/25000_att_4.graph"},
//    {
//    }
//  },
//  {
//    "grid_init",
//    BenchType::PROPORTION,
//    10,
//    1e-4,
//    0,
//    { "grid/Laplace_200x200.mtx" },
//    {
//      {"init", "time-init", "time-main", "none"},
//    }
//  },
//  {
//    "grid_stretch",
//    BenchType::STRETCH,
//    1,
//    1e-4,
//    0,
//    {"grid/Laplace_30x30.mtx"},
//    {}
//  },
//  {
//    "barabasi_stretch",
//    BenchType::STRETCH,
//    1,
//    1e-4,
//    0,
//    {"barabasi/1000_att_4.graph"},
//    {}
//  },
//  {
//    "grid_bench_2",
//    BenchType::ASYMPTOTIC,
//    3,
//    1e-4,
//    0,
//    {  "grid/Laplace_5x5.mtx",
//       "grid/Laplace_10x10.mtx",
//       "grid/Laplace_20x20.mtx",
//       "grid/Laplace_30x30.mtx",
//       "grid/Laplace_40x40.mtx",
//       "grid/Laplace_50x50.mtx",
//       "grid/Laplace_60x60.mtx",
//       "grid/Laplace_70x70.mtx",
//       "grid/Laplace_80x80.mtx",
//       "grid/Laplace_90x90.mtx",
//       "grid/Laplace_100x100.mtx",
//       "grid/Laplace_200x200.mtx",
//       "grid/Laplace_250x250.graph",
//       "grid/Laplace_300x300.mtx"},
//    {
//      {"time", "nodes", "time-mean", "time-var"},
//      {"cycles", "nodes", "cycles-mean", "cycles-var"},
//      {"flops", "nodes", "flops-mean", "flops-var"},
//      {"memory", "nodes", "memory-mean", "memory-var"},
//    }
//  },
//  {
//    "grid_bench_manual",
//    BenchType::ASYMPTOTIC,
//    10,
//    1e-4,
//    0,
//    GRIDS,
//    {
//      {"time", "nodes", "time-mean", "time-var"},
//      {"cycles", "nodes", "cycles-mean", "cycles-var"},
//      {"flops", "nodes", "flops-mean", "flops-var"},
//      {"memory", "nodes", "memory-mean", "memory-var"},
//    }
//  },
//  {
//    "grid3_bench",
//    BenchType::ASYMPTOTIC,
//    10,
//    1e-4,
//    0,
//    GRIDS_3D,
//    {
//      {"time", "nodes", "time-mean", "time-var"},
//      {"cycles", "nodes", "cycles-mean", "cycles-var"},
//      {"flops", "nodes", "flops-mean", "flops-var"},
//      {"memory", "nodes", "memory-mean", "memory-var"},
//    }
//  },
//  {
//    "barabasi_bench",
//    BenchType::ASYMPTOTIC,
//    10,
//    1e-4,
//    0,
//    BARABASI,
//    {
//      {"time", "nodes", "time-mean", "time-var"},
//      {"cycles", "nodes", "cycles-mean", "cycles-var"},
//      {"flops", "nodes", "flops-mean", "flops-var"},
//      {"memory", "nodes", "memory-mean", "memory-var"},
//    }
//  },
    {
      "manual",
      BenchType::ASYMPTOTIC,
      5,
      1e-4,
      0,
      {"../lesmis.graph", "../airfoil1.graph", "../PGPgiantcompo.graph", "luxembourg.osm.graph"},
      {
        {"time", "nodes", "time-mean", "time-var"},
        {"cycles", "nodes", "cycles-mean", "cycles-var"},
        {"flops", "nodes", "flops-mean", "flops-var"},
        {"memory", "nodes", "memory-mean", "memory-var"},
      }
    },
//  {
//    "grid_convergence",
//    BenchType::CONVERGENCE,
//    1,
//    1e-4,
//    0,
//    {"grid/Laplace_100x100.mtx"},
//    {{"residual_conv", "iters", "residual", "none"}, {"energy_conv", "iters", "energy", "none"}}
//  },
//  {
//    "grid_convergence_weighted",
//    BenchType::CONVERGENCE,
//    1,
//    1e-4,
//    0,
//    {"grid/Laplace_100x100_weighted.graph"},
//    {{"residual_conv", "iters", "residual", "none"}, {"energy_conv", "iters", "energy", "none"}}
//  },
//  {
//    "grid3_convergence",
//    BenchType::CONVERGENCE,
//    1,
//    1e-4,
//    0,
//    {"grid3/Laplace_20x20x20.graph"},
//    {
//      {"residual_conv", "iters", "residual", "none"},
//    }
//  },
//  {
//    "barabasi_convergence",
//    BenchType::CONVERGENCE,
//    1,
//    1e-4,
//    0,
//    {"barabasi/25000_att_4.graph"},
//    {{"residual_conv", "iters", "residual", "none"}, {"energy_conv", "iters", "energy", "none"}}
//  },
//  {
//    "barabasi_convergence_unweighted",
//    BenchType::CONVERGENCE,
//    1,
//    1e-4,
//    0,
//    {"barabasi/25000_att_4_unweighted.graph"},
//    {{"residual_conv", "iters", "residual", "none"}, {"energy_conv", "iters", "energy", "none"}}
//  },
//  {
//    "cong",
//    BenchType::CONVERGENCE,
//    1,
//    1e-2,
//    123456,
//    { "lesmis.graph" },
//    {}
//  },
};

// Available algorithms
using FlowAlgo = pair<string, function<Vector(const Graph&, const Vector&, SolverStatus&)>>;
vector<FlowAlgo> flow_algos = {
//    {"Uniform cycle, trivial flow, min. dist ST",    solveLaplacian<UniformCycleDistribution, TrivialFlow, minDistanceST>},
//    {"Uniform cycle, trivial flow, min. weight ST",  solveLaplacian<UniformCycleDistribution, TrivialFlow, minWeightST>},
//    {"Uniform cycle, trivial flow, low stretch ST", solveLaplacian<UniformCycleDistribution, TrivialFlow, lowStretchST>},
//    {"Uniform cycle, LCA flow, min. dist ST",        solveLaplacian<UniformCycleDistribution, LCAFlow, minDistanceST>},
//    {"Uniform cycle, LCA flow, min. weight ST",      solveLaplacian<UniformCycleDistribution, LCAFlow, minWeightST>},
//    {"Uniform cycle, LCA flow, low stretch ST",     solveLaplacian<UniformCycleDistribution, LCAFlow, lowStretchST>},
    {"Uniform cycle, log flow, min. dist ST",        solveLaplacian<UniformCycleDistribution, LogFlow, minDistanceST>},
    {"Uniform cycle, log flow, min. weight ST",      solveLaplacian<UniformCycleDistribution, LogFlow, minWeightST>},
//    {"Uniform cycle, log flow, special ST",      solveLaplacian<UniformCycleDistribution, LogFlow, specialGridST>},
//    {"Uniform cycle, log flow, low stretch ST",     solveLaplacian<UniformCycleDistribution, LogFlow, lowStretchST>},
//    {"Stretch cycle, trivial flow, min. dist ST",    solveLaplacian<StretchCycleDistribution, TrivialFlow, minDistanceST>},
//    {"Stretch cycle, trivial flow, min. weight ST",  solveLaplacian<StretchCycleDistribution, TrivialFlow, minWeightST>},
//    {"Stretch cycle, trivial flow, low stretch ST", solveLaplacian<StretchCycleDistribution, TrivialFlow, lowStretchST>},
//    {"Stretch cycle, LCA flow, min. dist ST",        solveLaplacian<StretchCycleDistribution, LCAFlow, minDistanceST>},
//    {"Stretch cycle, LCA flow, min. weight ST",      solveLaplacian<StretchCycleDistribution, LCAFlow, minWeightST>},
//    {"Stretch cycle, LCA flow, low stretch ST",     solveLaplacian<StretchCycleDistribution, LCAFlow, lowStretchST>},
    {"Stretch cycle, log flow, min. dist ST",        solveLaplacian<StretchCycleDistribution, LogFlow, minDistanceST>},
    {"Stretch cycle, log flow, min. weight ST",      solveLaplacian<StretchCycleDistribution, LogFlow, minWeightST>},
//    {"Stretch cycle, log flow, special ST",      solveLaplacian<StretchCycleDistribution, LogFlow, specialGridST>},
//    {"Stretch cycle, log flow, low stretch ST",     solveLaplacian<StretchCycleDistribution, LogFlow, lowStretchST>},
};

using CGAlgo = pair<string, pair<Vector,bool>(*)(const CompressedMatrix&, const Vector&, double)>;
vector<CGAlgo> cg_algos = {
//  {"CSR: CG, identity preconditioner", solveConjugateGradient<IdentityPreconditioner, CompressedMatrix>},
//  {"CSR: CG, Jacobi preconditioner", solveConjugateGradient<DiagonalPreconditioner, CompressedMatrix>},
};

#ifndef NOEIGEN
/* All combinations of eigen algorithms */
using EigenAlgo = pair<string, pair<EigenVector, int>(*)(const EigenMatrix&, const EigenVector&, edgeweight)>;
vector<EigenAlgo> eigen_algos = {
  {"Eigen: CG, identity preconditioner", solveEigen<EigenCG>},
  {"Eigen: CG, Jacobi preconditioner", solveEigen<EigenCGDiagPrecond>},
};
#endif

#ifndef NOPARALUTION
using ParalutionAlgo = pair<string, int(*)(const ParalutionMatrix&, const ParalutionVector&, ParalutionVector&, edgeweight)>;
vector<ParalutionAlgo> paralution_algos = {
  //{"Paralution: CG, identity preconditioner", solveParalution<ParalutionCG, ParalutionIdentity>}
  //{"Paralution: CG, Jacobi preconditioner", solveParalution<ParalutionCG, ParalutionJacobi>}
};
#endif

// Origins of graphs
const unordered_map<string, string> GRAPH_ORIGIN = {
  {"lesmis.graph", "http://www.cc.gatech.edu/dimacs10/archive/clustering.shtml"}
};

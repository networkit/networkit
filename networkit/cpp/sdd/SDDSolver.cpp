/*
 * SDDSolver.cpp
 *
 *  Created on: May 04, 2014
 *      Author: dhoske
 */

#include <vector>
#include <algorithm>
#include <atomic>

#include "Config.h"
#include "SDDSolver.h"

using namespace std;

namespace NetworKit {
namespace SDD {

bool isSDD(const Matrix& A) {
  if (!isSymmetric(A)) {
    return false;
  }

  /* Criterion: a_ii >= \sum_{j != i} a_ij */
  vector<double> row_sum(A.numberOfRows());
  A.parallelForNonZeroElementsInRowOrder([&] (node i, node j, double value) {
    if (i == j) {
      row_sum[i] += value;
    } else {
      row_sum[i] -= abs(value);
    }
  });
  return all_of(row_sum.begin(), row_sum.end(), [] (double val) {return val > -EPSILON;});
}

bool isLaplacian(const Matrix& A) {
  if (!isSymmetric(A)) {
    return false;
  }

  /* Criterion: \forall_i \sum_j A_ij = 0  */
  vector<double> row_sum(A.numberOfRows());
  std::atomic<bool> right_sign(true);
  A.parallelForNonZeroElementsInRowOrder([&] (node i, node j, double value) {
    if (i != j && value > EPSILON) {
      right_sign = false;
    }
    row_sum[i] += value;
  });
  return right_sign && all_of(row_sum.begin(), row_sum.end(), [] (double val) {return abs(val) < EPSILON;});
}

Graph laplacianToGraph(const Matrix& A) {
  Graph G(A.numberOfRows(), true);
  A.forNonZeroElementsInRowOrder([&] (node i, node j, double value) {
    if (i > j) {
      G.addEdge(i, j, -value);
    }
  });
//  for (index i = 0; i < A.numberOfRows(); ++i) {
//	  double sum = 0.0;
//	  double dEntry = 0.0;
//	  A.forNonZeroElementsInRow(i, [&](index i, index j, double value) {
//		  if (i > j) {
//			  G.addEdge(i, j, -value);
//			  sum += -value;
//		  } else if (i == j) {
//			  dEntry = value;
//		  }
//	  });
//
//	  if (dEntry != sum) { // self-loop
//		  G.addEdge(i, i, dEntry - sum);
//	  }
//  }




  return G;
}

Graph sddToLaplacian(const Matrix& A) {
  assert(isSDD(A));
  count n = A.numberOfColumns();

  /* Compute row sum and excess */
  vector<double> abs_row_sum(n);
  A.parallelForNonZeroElementsInRowOrder([&] (node i, node j, double value) {
    if (i != j) {
      abs_row_sum[i] += abs(value);
    }
  });

  vector<double> excess(n);
#pragma omp parallel for
  for (index i = 0; i < n; ++i) {
     excess[i] = A(i, i) - abs_row_sum[i];
  }

  /* Set lower and upper diagonal elements */
  Graph G(2*n, true);
  for (index i = 0; i < n; ++i) {
    G.addEdge(i, i + n, excess[i]/2);
  }

  /* Distribute original entries of A into the 2 times larger graph G */
  A.forNonZeroElementsInRowOrder([&] (node i, node j, double value) {
    if (i < j) {
      /* Off-diagonals */
      if (value < 0) {
        G.addEdge(i,     j,     -value);
        G.addEdge(i + n, j + n, -value);
      } else {
        G.addEdge(i + n,     j, value);
        G.addEdge(    i, j + n, value);
      }
    }
  });

  return G;
}

Vector sddToLaplacian(const Vector& b) {
  count n = b.getDimension();
  Vector y(2*n, 0);
  b.parallelForElements([&] (index i, double value) {
    y[i]     = value;
    y[i + n] = -value;
  });
  return y;
}

}
}

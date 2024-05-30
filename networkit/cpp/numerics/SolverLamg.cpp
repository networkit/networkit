/*
 * SolverLamg.cpp
 *
 *  Created on: 02.05.2024
 *      Author: Lukas
 */

#include <networkit/numerics/LAMG/SolverLamg.hpp>

namespace NetworKit {

template <>
void SolverLamg<CSRMatrix>::minRes(const index level, Vector &x, const Vector &r) const {
    if (numActiveIterates[level] > 0) {
        count n = numActiveIterates[level];

        std::vector<index> ARowIdx(r.getDimension() + 1);
        std::vector<index> ERowIdx(r.getDimension() + 1);

#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(r.getDimension()); ++i) {
            for (index k = 0; k < n; ++k) {
                double AEvalue = r[i] - rHistory[level][k][i];
                if (std::fabs(AEvalue) > 1e-25) {
                    ++ARowIdx[i + 1];
                }

                double Eval = history[level][k][i] - x[i];
                if (std::fabs(Eval) > 1e-25) {
                    ++ERowIdx[i + 1];
                }
            }
        }

        for (index i = 0; i < r.getDimension(); ++i) {
            ARowIdx[i + 1] += ARowIdx[i];
            ERowIdx[i + 1] += ERowIdx[i];
        }

        std::vector<index> AColumnIdx(ARowIdx[r.getDimension()]);
        std::vector<double> ANonZeros(ARowIdx[r.getDimension()]);

        std::vector<index> EColumnIdx(ERowIdx[r.getDimension()]);
        std::vector<double> ENonZeros(ERowIdx[r.getDimension()]);

#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(r.getDimension()); ++i) {
            for (index k = 0, aIdx = ARowIdx[i], eIdx = ERowIdx[i]; k < n; ++k) {
                double AEvalue = r[i] - rHistory[level][k][i];
                if (std::fabs(AEvalue) > 1e-25) {
                    AColumnIdx[aIdx] = k;
                    ANonZeros[aIdx] = AEvalue;
                    ++aIdx;
                }

                double Eval = history[level][k][i] - x[i];
                if (std::fabs(Eval) > 1e-25) {
                    EColumnIdx[eIdx] = k;
                    ENonZeros[eIdx] = Eval;
                    ++eIdx;
                }
            }
        }

        CSRMatrix AE(r.getDimension(), n, ARowIdx, AColumnIdx, ANonZeros, 0.0, true);
        CSRMatrix E(r.getDimension(), n, ERowIdx, EColumnIdx, ENonZeros, 0.0, true);

        Vector alpha = smoother.relax(CSRMatrix::mTmMultiply(AE, AE), CSRMatrix::mTvMultiply(AE, r),
                                      Vector(n, 0.0), 10);
        x += E * alpha;
    }
}

template <>
void SolverLamg<DynamicMatrix>::minRes(const index level, Vector &x, const Vector &r) const {
    if (numActiveIterates[level] > 0) {
        count n = numActiveIterates[level];

        DynamicMatrix AE(r.getDimension(), n);
        DynamicMatrix E(r.getDimension(), n);

#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(r.getDimension());
             ++i) {                                            // iterate rows of A and E: 0 to r
            for (index k = 0; k < n; ++k) {                    // iterate columns: 0 to n
                double AEvalue = r[i] - rHistory[level][k][i]; // value at (i,k)
                if (std::fabs(AEvalue)
                    > 1e-25) { // if there is a value, increase row count of row i+1 by 1
                    AE.setValue(i, k, AEvalue);
                }

                double Eval = history[level][k][i] - x[i];
                if (std::fabs(Eval) > 1e-25) {
                    E.setValue(i, k, Eval);
                }
            }
        }

        Vector alpha = smoother.relax(DynamicMatrix::mTmMultiply(AE, AE),
                                      DynamicMatrix::mTvMultiply(AE, r), Vector(n, 0.0), 10);
        x += E * alpha;
    }
}

template <>
void SolverLamg<DenseMatrix>::minRes(const index level, Vector &x, const Vector &r) const {
    if (numActiveIterates[level] > 0) {
        count n = numActiveIterates[level];

        DenseMatrix AE(r.getDimension(), n);
        DenseMatrix E(r.getDimension(), n);

#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(r.getDimension());
             ++i) {                                            // iterate rows of A and E: 0 to r
            for (index k = 0; k < n; ++k) {                    // iterate columns: 0 to n
                double AEvalue = r[i] - rHistory[level][k][i]; // value at (i,k)
                if (std::fabs(AEvalue)
                    > 1e-25) { // if there is a value, increase row count of row i+1 by 1
                    AE.setValue(i, k, AEvalue);
                }

                double Eval = history[level][k][i] - x[i];
                if (std::fabs(Eval) > 1e-25) {
                    E.setValue(i, k, Eval);
                }
            }
        }

        DenseMatrix AEtransposed = AE.transpose();

        Vector alpha = smoother.relax(AEtransposed * AE, AEtransposed * r, Vector(n, 0.0), 10);
        x += E * alpha;
    }
}

template class SolverLamg<CSRMatrix>;
template class SolverLamg<DenseMatrix>;
template class SolverLamg<DynamicMatrix>;

} // namespace NetworKit

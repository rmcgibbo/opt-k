#include <cmath>
#include <string>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <dlfcn.h>
#include "cycledeterminant.hpp"

extern "C" {
    // LAPACK LU decomposition, which we use to compute the determinant of our
    // matrix
    int dgetrf_(const int*, const int *, const double *, const int *, int *, int *);
}

namespace CycleDeterminant {


CycleMatrixBuilder::CycleMatrixBuilder(const int n_states)
    : n_states_(n_states),
      n_cycles_(0)
{
    for (int i = 1; i < n_states_; i++) {
        for (int j = i+1; j < n_states_; j++) {
            std::vector<int> cycle;
            cycle.push_back(indexInTriu(0, i)); // in "forward order" 0->i
            cycle.push_back(indexInTriu(i, j)); // in "forward order" i->j
             // BACKWARD! this edge actually corresponds to j->0
            cycle.push_back(indexInTriu(0, j));
            cycles_.push_back(cycle);
        }
    }
    n_cycles_ = cycles_.size();
}

CycleMatrixBuilder::~CycleMatrixBuilder() {
}

double CycleMatrixBuilder::logSqrtDetCycleMatrix(const double* x) {
    std::vector<double> A(n_cycles_ * n_cycles_);
    for (int i = 0; i < n_cycles_; i++) {
        for (int j = 0; j < n_cycles_; j++) {
            if (i == j) {
                A[i*n_cycles_ + i] = \
                    exp(-x[cycles_[i][0]]) +  exp(-x[cycles_[i][1]]) +  exp(-x[cycles_[i][2]]);
            } else {
                // Look for edges that are in the intersection set of the two
                // cycles
                
                // the 0th and 1st edges in both cycle are forward order, so
                // if they intersect the sign of the term will be positive
                for (int a = 0; a < 2; a++)
                    for (int b = 0; b < 2; b++)
                        if (cycles_[i][a] == cycles_[j][b])
                            A[i*n_cycles_ + j] += exp(-x[cycles_[i][a]]);

                // The 2nd edge is in backward order, so if two 2nd edges
                // intersect, they will also have the same sign.
                if (cycles_[i][2] == cycles_[j][2])
                    A[i*n_cycles_ + j] += exp(-x[cycles_[i][2]]);

                // Since the 0th and 1st are in forward order but the 2nd is
                // backward, if there's an intersection, they will have opposite
                // directions and get a negative sign.
                for (int a = 0; a < 2; a++) {
                    if (cycles_[i][a] == cycles_[j][2])
                        A[i*n_cycles_ + j] -= exp(-x[cycles_[i][a]]);
                    if (cycles_[i][2] == cycles_[j][a])
                        A[i*n_cycles_ + j] -= exp(-x[cycles_[j][a]]);
                }
            }
        }
    }

    double value = 0.5*logdet(A);
    if (value != value)
        throw std::overflow_error("NaN determinant of A(x) in cycledeterminant.cpp");
    return value;
}


double CycleMatrixBuilder::logdet(std::vector<double> const& A) {
    int N = round(sqrt(A.size()));
    if (N*N != A.size())
        throw std::runtime_error("A must be N by N");

    int info = 0;
    std::vector<int> pivots(N);
    dgetrf_(&N, &N, &A[0], &N, &pivots[0], &info);

    int neg = 0;
    double logdet = 0.0;
    for (int c1 = 0; c1 < N; c1++) {
        logdet += log(std::abs(A[N*c1 + c1]));
        if (A[N*c1 + c1] < 0) neg = !neg;
        if (pivots[c1] != (c1+1)) neg = !neg;
    }
    
    // the determinant of the cycle matrix is NOT supposed to be negative
    if (neg)
        throw std::runtime_error("Negative determiantn!!");

    return logdet;
}

int CycleMatrixBuilder::indexInTriu(int i, int j) {
    // Compte the index in the 1D symmetric storage for the edge from i to j
    if (i > j)
        throw std::runtime_error("Error");

    return static_cast<int>(-0.5 * i*i + (n_states_ - 0.5) * i + j);
}

// close namespace
};


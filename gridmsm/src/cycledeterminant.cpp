#include <cmath>
#include <string>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <dlfcn.h>
#include "cycledeterminant.hpp"
namespace CycleDeterminant {


CycleMatrixBuilder::CycleMatrixBuilder(const int n_states, std::string lapack_lite_lib)
    : n_states_(n_states),
      n_cycles_(0)
{
    for (int i = 1; i < n_states_; i++) {
        for (int j = i+1; j < n_states_; j++) {
            std::vector<int> cycle;
            cycle.push_back(indexInTriu(0, i));
            cycle.push_back(indexInTriu(i, j));
            cycle.push_back(indexInTriu(0, j));
            cycles_.push_back(cycle);
        }
    }
    n_cycles_ = cycles_.size();

    //void* pyhandle = dlopen(python_lib.c_str(),  RTLD_LAZY | RTLD_LOCAL);
    //if (!pyhandle)
    //    throw std::runtime_error(dlerror());
    lapack_lite = dlopen(lapack_lite_lib.c_str(), RTLD_LAZY | RTLD_LOCAL);
    if (!lapack_lite)
        throw std::runtime_error(dlerror());
    //dlclose(pyhandle);
    dgetrf_ = (int (*)(const int*, const int *, const double *, const int *, int *, int *)) dlsym(lapack_lite, "dgetrf_");
    if (!dgetrf_)
        throw std::runtime_error(dlerror());
}

CycleMatrixBuilder::~CycleMatrixBuilder() {
    dlclose(lapack_lite);
}

double CycleMatrixBuilder::logSqrtDetCycleMatrix(const double* u) {
    std::vector<double> A(n_cycles_ * n_cycles_);
    for (int i = 0; i < n_cycles_; i++) {
        for (int j = 0; j < n_cycles_; j++) {
            if (i == j) {
                A[i*n_cycles_ + i] = \
                    exp(-u[cycles_[i][0]]) +  exp(-u[cycles_[i][1]]) +  exp(-u[cycles_[i][2]]);
            } else {
                for (int x = 0; x < 2; x++)
                    for (int y = 0; y < 2; y++)
                        if (cycles_[i][x] == cycles_[j][y])
                            A[i*n_cycles_ + j] += exp(-u[cycles_[i][x]]);
                for (int x = 0; x < 2; x++) {
                    if (cycles_[i][x] == cycles_[j][2])
                        A[i*n_cycles_ + j] -= exp(-u[cycles_[i][x]]);
                    if (cycles_[i][2] == cycles_[j][x])
                        A[i*n_cycles_ + j] -= exp(-u[cycles_[j][x]]);
                }
            }
        }
    }
    //
    // cout << "[ ";
    // for (int i = 0; i < n_cycles_; i++) {
    //     cout << "[ ";
    //     for (int j = 0; j < n_cycles_; j++)
    //         std::cout << A[i*n_cycles_ + j] << ", ";
    //     std::cout << "]," << std::endl;
    // }
    // cout << "]";

    return 0.5*logdet(A);
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

    return logdet;
}

int CycleMatrixBuilder::indexInTriu(int i, int j) {
    if (i > j)
        throw std::runtime_error("Error");

    return static_cast<int>(-0.5 * i*i + (n_states_ - 0.5) * i + j);
}

// close namespace
};


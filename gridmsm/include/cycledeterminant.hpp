#include <vector>

namespace CycleDeterminant {
class CycleMatrixBuilder {

public:
    CycleMatrixBuilder(const int n_states, std::string lapack_lite_lib);
    double logSqrtDetCycleMatrix(const double* u);
    ~CycleMatrixBuilder();

private:
    int indexInTriu(int i, int j);
    double logdet(std::vector<double> const& A);

    const int n_states_;
    int n_cycles_;
    void* lapack_lite;
    std::vector<std::vector<int> > cycles_;
    int (*dgetrf_)(const int*, const int *, const double *, const int *, int *, int *);
};

}; // closes namespace
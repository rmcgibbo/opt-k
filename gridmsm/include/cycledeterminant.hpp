#include <vector>

namespace CycleDeterminant {
class CycleMatrixBuilder {

public:
    CycleMatrixBuilder(const int n_states);
    double logSqrtDetCycleMatrix(const double* u);
    ~CycleMatrixBuilder();

private:
    int indexInTriu(int i, int j);
    double logdet(std::vector<double> const& A);

    const int n_states_;
    int n_cycles_;
    std::vector<std::vector<int> > cycles_;
};

}; // closes namespace
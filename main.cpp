#include <iostream>
#include <vector>
#include <cassert>
#include <iomanip>
#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"
#include "CoinMpsIO.hpp"
#include "slu_ddefs.h"


class InnerPointSolver{
public:
    explicit InnerPointSolver();
    ~InnerPointSolver();

    bool loadMPSMatrix(const std::string &filename);
    void solve();

private:
    void printCoinPackedMatrix(const CoinPackedMatrix& T);
    void printSuperMatrix(const SuperMatrix &A);
    CoinPackedMatrix createDiagonalMatrix(const std::vector<double>& diagonalValues);
    CoinPackedMatrix multiplyMatrices(const CoinPackedMatrix& A, const CoinPackedMatrix& B);
    void convertToSuperLU(CoinPackedMatrix& data, SuperMatrix& A_slu);
    double calculateLambda(const std::vector<double>& x, const std::vector<double>& s, const std::vector<double>& r, double gamma);


    CoinPackedMatrix A;
    size_t Arows, Acols;
    std::vector<double> c, b, x, r, g;
    const double gamma = 0.3333;
    double lambda;

    CoinPackedMatrix D, A_T, AD, ADA_t;
    SuperMatrix T_super, b_super, L, U;
    int *perm_c, *perm_r;
    superlu_options_t options;
    SuperLUStat_t stat;

    std::vector<double> Ax, d, ADc, A_Tu, q, u, s;
    double res;
    bool solved;
};

InnerPointSolver::InnerPointSolver() {
    Arows = 0;
    Acols = 0;
    lambda = 0;
    solved = false;
}

InnerPointSolver::~InnerPointSolver() {
    Destroy_SuperMatrix_Store(&T_super);
    Destroy_SuperMatrix_Store(&b_super);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    delete[] perm_c;
    delete[] perm_r;
}


void InnerPointSolver::printCoinPackedMatrix(const CoinPackedMatrix& T) {
    int numCols = T.getNumCols();
    const int* colStarts = T.getVectorStarts();
    const int* rowIndices = T.getIndices();
    const double* values = T.getElements();

    std::cout << "Matrix T:" << std::endl;
    for (int col = 0; col < numCols; ++col) {
        int start = colStarts[col];
        int end = colStarts[col + 1];

        for (int idx = start; idx < end; ++idx) {
            int row = rowIndices[idx];
            double value = values[idx];
            std::cout << "T(" << row << ", " << col << ") = " << value << std::endl;
        }
    }
}


void InnerPointSolver::printSuperMatrix(const SuperMatrix &A) {
    if (A.Stype != SLU_NC) {
        std::cerr << "Matrix is not in form if Compressed Sparse Column (CSC)." << std::endl;
        return;
    }
    if (A.Dtype != SLU_D) {
        std::cerr << "Matrix is not of double." << std::endl;
        return;
    }
    if (A.Mtype != SLU_GE) {
        std::cerr << "Matrix is not general." << std::endl;
        return;
    }

    NCformat *store = (NCformat *)A.Store;
    double *values = (double *)store->nzval;
    int *rowIndex = store->rowind;
    int *colPtr = store->colptr;

    int numRows = A.nrow;
    int numCols = A.ncol;

    std::cout << "SuperMatrix (CSC Format):" << std::endl;
    std::cout << "Rows: " << numRows << ", Columns: " << numCols << std::endl;

    for (int col = 0; col < numCols; ++col) {
        std::cout << "Column " << col << ":" << std::endl;
        for (int idx = colPtr[col]; idx < colPtr[col + 1]; ++idx) {
            int row = rowIndex[idx];
            double value = values[idx];
            std::cout << "  A(" << row << ", " << col << ") = " << value << std::endl;
        }
    }
}


CoinPackedMatrix InnerPointSolver::createDiagonalMatrix(const std::vector<double>& diagonalValues) {
    size_t size = diagonalValues.size();

    CoinPackedMatrix matrix;

    for (size_t i = 0; i < size; ++i) {
        if (diagonalValues[i] != 0.0) {
            int index[] = {i};
            double value[] = {diagonalValues[i]};
            matrix.appendCol(1, index, value);
        } else {
            int index[] = {};
            double value[] = {};
            matrix.appendCol(0, index, value);
        }
    }

    return matrix;
}


CoinPackedMatrix InnerPointSolver::multiplyMatrices(const CoinPackedMatrix& A, const CoinPackedMatrix& B) {
    if (A.getNumCols() != B.getNumRows()) {
        std::cerr << "Number of columns of A must be equal number of rows of B!" << std::endl;
        exit(1);
    }

    int m = A.getNumRows();
    int p = B.getNumCols();
    CoinPackedVector rowA, colB;
    CoinPackedMatrix C;
    CoinPackedMatrix A_temp = A;
    A_temp.reverseOrdering();

    for (int j = 0; j < p; ++j) {
        colB = B.getVector(j);
        int* indicesB = colB.getIndices();
        double* valuesB = colB.getElements();
        int lengthB = colB.getNumElements();

        std::vector<int> rowIndices;
        std::vector<double> rowValues;

        for (int i = 0; i < m; ++i) {
            double dotProduct = 0.0;

            rowA = A_temp.getVector(i);
            int* indicesA = rowA.getIndices();
            double* valuesA = rowA.getElements();
            int lengthA = rowA.getNumElements();

            for (int k = 0; k < lengthA; ++k) {
                for (int l = 0; l < lengthB; ++l) {
                    if (indicesA[k] == indicesB[l]) {
                        dotProduct += valuesA[k] * valuesB[l];
                    }
                }
            }

            if (dotProduct != 0.0) {
                rowIndices.push_back(i);
                rowValues.push_back(dotProduct);
            }
        }

        if (!rowIndices.empty()) {
            C.appendCol(rowIndices.size(), rowIndices.data(), rowValues.data());
        } else {
            C.appendCol(0, nullptr, nullptr);
        }
    }

    return C;
}


bool InnerPointSolver::loadMPSMatrix(const std::string &filename) {
    CoinMpsIO mpsReader;
    int numErr = mpsReader.readMps(filename.c_str(), "mps");
    if (numErr != 0) {
        std::cerr << "Error loading file " << filename << std::endl;
        return false;
    }

    A = *mpsReader.getMatrixByCol();
    const double *cArray = mpsReader.getObjCoefficients();
    const double *bArray = mpsReader.getRightHandSide();

    c = std::vector<double>(cArray, cArray + A.getNumCols());
    b = std::vector<double>(bArray, bArray + A.getNumRows());

    Acols = A.getNumCols();
    Arows = A.getNumRows();

    std::cout << "Matrix A loaded. Rows: " << A.getNumRows()
              << ", Columns: " << A.getNumCols() << std::endl;


    return true;
}

void InnerPointSolver::convertToSuperLU(CoinPackedMatrix& data, SuperMatrix& A_slu) {
    int numRows = data.getNumRows();
    int numCols = data.getNumCols();
    const double* values = data.getElements();
    const int* rowIndices = data.getIndices();
    const int* colStarts = data.getVectorStarts();

    A_slu.Store = new NCformat;
    NCformat* store = (NCformat*)A_slu.Store;

    store->nzval = const_cast<double*>(values);
    store->rowind = const_cast<int*>(rowIndices);
    store->colptr = const_cast<int*>(colStarts);

    store->nnz = colStarts[numCols];
    A_slu.nrow = numRows;
    A_slu.ncol = numCols;

    A_slu.Stype = SLU_NC;
    A_slu.Dtype = SLU_D;
    A_slu.Mtype = SLU_GE;
}

void cleanupSuperLU(SuperMatrix& A_slu) {
    Destroy_SuperMatrix_Store(&A_slu);
}

double InnerPointSolver::calculateLambda(const std::vector<double>& x, const std::vector<double>& s, const std::vector<double>& r, double gamma) {
    double minRatio = 1e10;
    double ratio = 0.;
    double r_length = 0.;

    for (size_t j = 0; j < x.size(); ++j) {
        if (s[j] < 0) {
            ratio = -x[j] / s[j];
            if (ratio < minRatio) {
                minRatio = ratio;
            }
        }
    }

    for (size_t i = 0; i < r.size(); ++i) {
        r_length += r[i] * r[i];
    }

    std::cout << "(r_length)^2 = " << r_length << std::endl;

    double lambdaBar = gamma * minRatio;

    return (r_length < 0.001) ? lambdaBar : std::min(1.0, lambdaBar);
}


void InnerPointSolver::solve() {
    x = std::vector<double>(Acols, 1.);
    r = std::vector<double>(Arows, 0.);
    g = std::vector<double>(Acols, 0.);
    Ax = std::vector<double>(Arows, 0.);
    d = std::vector<double>(Acols, 1.0);
    ADc = std::vector<double>(Arows, 0.);
    A_Tu = std::vector<double>(Acols, 0.);
    q = std::vector<double>(Arows, 1.);
    u = std::vector<double>(Arows, 0.);
    s = std::vector<double>(Acols, 0.);


    A.times(x.data(), Ax.data());

    for (size_t i = 0; i < A.getNumRows(); ++i) {
        r[i] = b[i] - Ax[i];
    }


    std::cout << "b = ";
    for (size_t i = 0; i < b.size(); ++i) {
        std::cout << b[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "c = ";
    for (size_t i = 0; i < c.size(); ++i) {
        std::cout << c[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "r = ";
    for (size_t i = 0; i < r.size(); ++i) {
        std::cout << r[i] << " ";
    }
    std::cout << std::endl;




    set_default_options(&options);
    StatInit(&stat);


    A_T = A;
    A_T.transpose();
    A_T.reverseOrdering();

    std::cout << "A:";
    printCoinPackedMatrix(A);
    std::cout << "A_T:";
    printCoinPackedMatrix(A_T);

    while (!solved) {
        for (size_t i = 0; i < d.size(); ++i) {
            d[i] = x[i] * x[i];
        }

        D = createDiagonalMatrix(d);

        std::cout << "D:";
        printCoinPackedMatrix(D);


        AD = multiplyMatrices(A, D);
        std::cout << "AD:";
        printCoinPackedMatrix(AD);

        ADA_t = multiplyMatrices(AD, A_T);

        std::cout << "ADA_t:";
        printCoinPackedMatrix(ADA_t);

        convertToSuperLU(ADA_t, T_super);

        std::cout << "T_super:";
        printSuperMatrix(T_super);

        AD.times(c.data(), ADc.data());


        for (size_t i = 0; i < q.size(); ++i) {
            q[i] = r[i] + ADc[i];
        }


        dCreate_Dense_Matrix(&b_super, Arows, 1, q.data(), Arows, SLU_DN, SLU_D, SLU_GE);

        perm_c = new int[Arows];
        perm_r = new int[Arows];

        int info;
        dgssv(&options, &T_super, perm_c, perm_r, &L, &U, &b_super, &stat, &info);

        if (info == 0) {
            double *sol = (double *)((DNformat *)b_super.Store)->nzval;
            std::copy(sol, sol + Arows, u.begin());

            std::cout << "Solution u:" << std::endl;
            for (double val : u) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        } else {
            std::cerr << "Error solving system, info: " << info << std::endl;
        }


        A_T.times(u.data(), A_Tu.data());

        for (size_t i = 0; i < g.size(); ++i) {
            g[i] = c[i] - A_Tu[i];
        }

        D.times(g.data(), s.data());

        for (size_t i = 0; i < s.size(); ++i) {
            s[i] = -s[i];
        }



        lambda = calculateLambda(x, s, r, gamma);

        std::cout << "Lambda: " << lambda << std::endl;

        for (size_t i = 0; i < x.size(); ++i) {
            x[i] += lambda * s[i];
        }

        for (size_t i = 0; i < r.size(); ++i) {
            r[i] *= (1 - lambda);
        }


        res = 0.;
        for (size_t i = 0; i < c.size(); ++i) {
            res += c[i]*x[i];
        }

        std::cout << "Result c*x = " << std::fixed << std::setprecision(6) << res << std::endl;

        double r_length = 0.;
        for (size_t i = 0; i < r.size(); ++i) {
            r_length += r[i] * r[i];
        }

        if (r_length < (1e-2)) {
            solved = true;
        }
    }

    std::cout << "Result x = ";
    for (size_t i = 0; i < x.size(); ++i) {
        std::cout << x[i] << " ";
    }
    std::cout << '\n';
}


int main() {
    InnerPointSolver IPS;
    if (!IPS.loadMPSMatrix("test.mps")) {
        return 1;
    }
    IPS.solve();

    return 0;
}

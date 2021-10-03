#pragma once

//
// Numerical analysis library v0.0.1
//

#include <iomanip>
#include <iostream>

#include <stdexcept>

#include <iterator>
#include <vector>

#include <algorithm>
#include <random>

#include <cmath>

#include <initializer_list>

double real_rand();
int int_rand();

double real_rand() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dis(0, 1.0);
    return dis(gen);
}

int int_rand() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_int_distribution<> dis(0, 9);
    return dis(gen);
}

template <class T = double>
class matrix {
    std::vector<T> data;
    int cols, rows;

public:
    matrix(int h, int w) : data(w * h, 0) {
        cols = w;
        rows = h;
    }

    matrix(int h, int w, T value) : matrix(w, h) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                this->set(i, j) = value;
            }
        }
    }

    matrix(std::initializer_list<std::initializer_list<T>> l)
        : matrix(l.size(), l.begin()->size()) {
        int count = 0;
        for (auto i : l) {
            for (auto j : i) {
                data[count++] = j;
            }
        }
    }

public:
    auto &get_data() { return data; }

    const auto &get_data() const { return data; }
    const auto &get(int y, int x) const { return data.at(x + y * cols); }

    const auto get_rows() const { return rows; }
    const auto get_cols() const { return cols; }

    auto &set(int y, int x, const T &val) { return data.at(x + y * cols) = val; }
    auto &set(int y, int x) { return data.at(x + y * cols); }

    friend std::ostream &operator<<(std::ostream &str, const matrix &m) {
        for (size_t i = 0; i < m.rows; i++) {
            std::string delimiter = "";
            str << "| " << std::showpoint << std::setprecision(6);
            for (int j = 0; j < m.cols; ++j) {
                str << delimiter << std::setw(12) << m.data.at(j + i * m.cols);
                delimiter = " , ";
            }
            str << " |\n";
        }
        return str;
    }

    friend matrix<T> operator*(const matrix &A, const matrix &B) {
        return dot_product(A, B);
    }

    friend matrix<T> operator+(matrix A, const matrix B) {
        if((A.get_cols() != B.get_cols()) || (A.get_rows() != B.get_rows())) {
            throw std::invalid_argument("Wrong matrix dimensions.");
        } else {
            auto &v1 = A.get_data();
            auto &v2 = B.get_data();
            for(int i = 0; i < v1.size(); i++)
                v1[i] += v2[i];

            return A;
        }
    }

    friend matrix<T> operator-(matrix A, const matrix B) {
        if((A.get_cols() != B.get_cols()) || (A.get_rows() != B.get_rows())) {
            throw std::invalid_argument("Wrong matrix dimensions.");
        } else {
            auto &v1 = A.get_data();
            auto &v2 = B.get_data();
            for(int i = 0; i < v1.size(); i++) {
                v1[i] -= v2[i];
            }
            return A;
        }
    }
};

// Generate matrix for storing roots.
template <class T = double>
auto X_matrix(const matrix<T> &S) {
    return matrix<T>(S.get_rows(), 1);
}

matrix<double> double_rand(int h, int w) {
    matrix<double> m(h, w);
    auto &data = m.get_data();
    std::generate(data.begin(), data.end(), real_rand);
    return m;
}

template <class T = double>
auto dot_product(const matrix<T> &A, const matrix<T> &B) {
    if (A.get_cols() != B.get_rows()) {
        throw std::invalid_argument("Wrong matrix dimensions.");
    }

    auto C = matrix(A.get_rows(), B.get_cols());

    for (int i = 0; i < C.get_rows(); ++i) {
        for (int j = 0; j < C.get_cols(); ++j) {
            for (int l = 0; l < A.get_cols(); ++l) {
                C.set(i, j) += (A.get(i, l) * B.get(l, j));
            }
        }
    }

    return C;
}

// Stack horizontally
template <class T = double>
matrix<T> hstack(const matrix<T> A, const matrix<T> B) {
    if (A.get_rows() != B.get_rows()) {
        throw std::invalid_argument("Rows are not equal.");
    }

    matrix<T> C(A.get_rows(), A.get_cols() + B.get_cols());

    // Stack first matrix
    for (int i = 0; i < A.get_rows(); ++i) {
        for (int j = 0; j < A.get_cols(); ++j) {
            C.set(i, j) = A.get(i, j);
        }
    }

    // Stack second matrix
    for (int i = 0; i < B.get_rows(); ++i) {
        for (int j = 0; j < B.get_cols(); ++j) {
            C.set(i, A.get_cols() + j) = B.get(i, j);
        }
    }

    return C;
}

template <class T = double>
auto transpose(const matrix<T> &M) {
    const int &M_rows = M.get_cols();
    const int &M_cols = M.get_cols();
    matrix<T> tM(M_cols, M_rows);
    for (int i = 0; i < M_cols; ++i) {
        for (int j = 0; j < M_rows; ++j) {
            tM.set(i, j) = M.get(j, i);
        }
    }
    return tM;
}

// Transform matrix into upper triangular form for gauss method
// M = (A|B) system of equations
// is_verbose - flag to output every step of transformation
template <class T = double>
auto to_triangle(matrix<T> M, bool is_verbose = false) {
    const int &rows = M.get_rows();
    const int &cols = M.get_cols();

    for (int k = 0; k < rows - 1; k++) {
        double n1 = M.get(k, k);
        for (int i = 1 + k; i < rows; i++) {
            double n2 = M.get(i, k);
            for (int j = 0; j < cols; j++) {
                M.set(i, j) = n1 * M.get(i, j) - n2 * M.get(k, j);
            }
        }
        if(is_verbose) std::cout << (k + 1) << ":\n" << M << "\n";
    }
    return M;
}

template <class T = double>
auto back_substitution(matrix<T> M) {
    const int &rows = M.get_rows();
    const int &cols = M.get_cols();

    matrix<T> X(rows, 1);

    for (int i = 0; i < rows; i++) {
        for (int j = 1; j <= i; ++j) {
            M.set((rows - 1) - i, cols - 1) -= (M.get((rows - 1) - i, (rows) - j) * X.get((rows) - j, 0));
        }
        X.set((rows - 1) - i, 0) = M.get((rows - 1) - i, cols - 1) / M.get((rows - 1) - i, (rows - 1) - i);
    }

    return X;
}

template <class T = double>
auto forward_substitution(matrix<T> M) {
    const int &rows = M.get_rows();
    const int &cols = M.get_cols();

    matrix<T> X(rows, 1);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < i; ++j) {
            M.set(i, cols - 1) -= (M.get(i, j) * X.get(j, 0));
        }
        X.set(i, 0) = M.get(i, cols - 1) / M.get(i, i);
    }

    return X;
}

template<class T = double>
bool is_symmetrical(const matrix<T> &M) {
    const int &rows = M.get_rows();
    for (int i = 0; i < rows; ++i) {
        for (int j = i + 1; j < rows; ++j) {
            if(M.get(i, j) != M.get(j, i))
                return false;
        }
    }
    return true;
}

// Return upper triangle matrix (including diagonals)
template<class T = double>
auto upper_triangle(matrix<T> M) {
    const int &rows = M.get_rows();
    for (int i = 0; i < rows; ++i) {
        for (int j = i + 1; j < rows; ++j) {
            M.set(j, i) = 0;
        }
    }
    return M;
}

template<class T = double>
auto lower_triangle(matrix<T> M) {
    return transpose(upper_triangle(M));
}


// Cholesky decomposition
template<class T = double>
auto cholesky(matrix<T> A) {
    matrix<T> L(A.get_rows(), A.get_cols(), 0);

    for (int i = 0; i < L.get_rows(); i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0;
            for (int k = 0; k < j; k++)
                sum += L.get(i, k) * L.get(j, k);
            if (i == j)
                L.set(i, j) = std::sqrt(A.get(i, i) - sum);
            else
                L.set(i, j) = (1.0 / L.get(j, j) * (A.get(i, j) - sum));
        }
    }

    return L;
}


// Solve system of linear equations using gauss method
template <class T = double>
auto gauss(matrix<T> M, bool is_verbose = false) {
    return back_substitution(to_triangle(M, is_verbose));
}

// Solve system of linear equations using square roots method 
template<class T = double>
matrix<T> sqrt_method(matrix<T> A, matrix<T> B) {

    matrix<T> L = cholesky(A);

    std::cout << "L:\n" << L << "\n";

    // T'y = B
    matrix<T> Y = forward_substitution(hstack(L, B));

    std::cout << "Y:\n" << Y << "\n";

    // Tx=y
    matrix<T> X = back_substitution(hstack(transpose(L), Y));

    return X;
}
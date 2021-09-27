#pragma once

//
// Numerical analysis library v0.0.1
//

#include <iostream>
#include <iomanip>

#include <stdexcept>

#include <vector>
#include <iterator>

#include <algorithm>
#include <random>

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



template<class T = double>
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
        : matrix( l.size(), l.begin()->size()) {
        int count = 0;
        for(auto i : l) {
            for(auto j : i) {
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
        for(size_t i = 0; i < m.rows; i++) {
            std::string delimiter = "";
            str << "| "
                << std::fixed << std::setprecision(5);
            for (int j = 0; j < m.cols; ++j) {
                str << delimiter
                    << std::setw(10)
                    << m.data.at(j + i * m.cols);
                delimiter = " , ";
            }
            str << " |\n";
        }
        return str;
    }

    void to_triangle() {

        auto &M = *this;

        for(int k = 0; k < rows - 1; k++) {
            double n1 = M.get(k, k);
            for(int i = 1 + k; i < rows; i++) {
                double n2 = M.get(i, k);
                for(int j = 0; j < cols; j++) {
                    M.set(i, j) = n1 * M.get(i, j) - n2 * M.get(k, j);
                }
            }
        }
    }

    friend matrix operator*(const matrix &A, const matrix &B) {
        return dot_product(A, B);
    }

};

matrix<double> randMatrix(int h, int w) {
    matrix<double> m(h, w);
    auto &data = m.get_data();
    std::generate(data.begin(), data.end(), int_rand);
    return m;
    // return std::move(m); // ?
}

template<class T = double>
matrix<T> dot_product(const matrix<T> &A, const matrix<T> &B) {
    if(A.get_cols() != B.get_rows()) {
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
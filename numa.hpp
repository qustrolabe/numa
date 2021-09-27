#pragma once

//
// Numerical analysis library v0.0.1
//

#include <iostream>
#include <iomanip>

#include <vector>
#include <iterator>

#include <algorithm>
#include <random>



double real_rand() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dis(-1.0, 1.0);
    // static std::uniform_int_distribution<> dis(0, 3);
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

    matrix(int w, int h, T value) : matrix(w, h) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                this->set(i, j) = value;
            }
        }
    }

public:
    auto &get_data() { return data; }
    const auto &get_data() const { return data; }
    const auto &get(int y, int x) const { return data.at(x + y * cols); }
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

    // DOT PRODUCT OF TWO MATRICES
    // TODO: MOVE TO EXTERNAL FUNCTION
    friend matrix operator*(const matrix &A, const matrix &B) {
        // A.cols == B.rows
        // C {A.rows, B.cols}

        if(A.cols != B.rows) {
            //Throw error
            std::cout << "Wrong demensions!!\n";
        }

        auto C = matrix(A.rows, B.cols);

        for (int i = 0; i < C.rows; ++i) {
            for (int j = 0; j < C.cols; ++j) {
                for (int l = 0; l < A.cols; ++l) {
                    C.set(i, j) += (A.get(i, l) * B.get(l, j));
                }
            }
        }

        std::cout << C << "\n";

        return matrix(1, 1);
    }

};

matrix<double> randMatrix(int h, int w) {
    matrix<double> m(h, w);
    auto &data = m.get_data();
    std::generate(data.begin(), data.end(), real_rand);
    return m;
    // return std::move(m); // ?
}
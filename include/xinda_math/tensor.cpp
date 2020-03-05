//
// Created by sunduda on 3/4/20.
//

#include "tensor.h"
#include <iostream>
#include <algorithm>

using namespace xinda_math;

template<std::size_t D, typename T>
template<typename cT>
tensor<D, T>::tensor(cT (&sizes)[D], T &value) {
    if (!(std::is_same<cT, std::size_t>::value || std::is_same<cT, unsigned short>::value ||
          std::is_same<cT, unsigned>::value || std::is_same<cT, unsigned long>::value))
        throw std::invalid_argument("Invalid typename!");
    reset(sizes, value);
}

template<std::size_t D, typename T>
template<typename cT>
void tensor<D, T>::reset(cT (&sizes)[D], T &value) {
    if (!(std::is_same<cT, std::size_t>::value || std::is_same<cT, unsigned short>::value ||
          std::is_same<cT, unsigned>::value || std::is_same<cT, unsigned long>::value))
        throw std::invalid_argument("Invalid typename!");
    cT(&new_sizes)[D - 1] = sizes + 1;
    tensor<D - 1, T> &temp_tensor(new_sizes, value);
    for (cT i = 0; i < sizes[0]; i++) this->mTensor.push_back(temp_tensor.get_value());
}

/*
template <typename T>
struct multidimensional_vector<0, T>
{
    typedef T type;
};

// Constructor
template <size_t D, typename T>
tensor<D, T>::tensor(size_t (&size)[D], T &value) {
}

template <size_t D, typename T>
bool tensor<D, T>::reset(size_t(&size)[D]) {
    size_t dims = 0;

    for (size_t i = 0; i < D; i++) {
        if (size[i] == 0) size[i] = 1;
    }

    try {
        while (true) {
            if (dims == D) break;
            if (dims == D - 1)
                mTensor.resize(size[dims], (T)0);
            else
                mTensor.resize(size[dims], multidimensional_vector<D - 1 - dims, T>({ 0 }));
            main_matrix_.resize(n_rows_, std::vector<T>(n_cols_, (T)0));
            main_matrix_.shrink_to_fit();
            // The size of compressed LU matrix is the same as its original matrix.
            matrix_lu_.resize(n_rows_, std::vector<double>(n_cols_, (double)0));
            matrix_lu_.shrink_to_fit();
            // Elementary Row/Column Operations mapping arrays initialisation.
            matrix_pr_.resize(n_rows_, (int)0);
            matrix_pc_.resize(n_cols_, (int)0);
            matrix_pr_.shrink_to_fit();
            matrix_pc_.shrink_to_fit();
        }
        catch (std::bad_alloc const&) {
                std::cout << "Memory allocation failed!" << std::endl;
                return false;
        }
        for (long i = 0; i < n_rows_; i++) matrix_pr_[i] = i;
        for (long j = 0; j < n_cols_; j++) matrix_pc_[j] = j;

        display_width_ = 4;
        FlagReset();

        return true;
    }
}

template <size_t D, typename T>
void tensor<D, T>::FlagReset() {
    // Set flags to initial state.
    is_decomposed_ = false;
    is_invertible_ = true;
}

template <size_t D, typename T>
void tensor<D, T>::SetDisplayWidth(const int& pcs) {
    int n_intd = 1;
    int max_n_intd = 1;
    for (int i = 0; i < main_matrix_.size(); i++) {
        for (int j = 0; j < main_matrix_[i].size(); j++) {
            n_intd = (unsigned)(log10(abs(main_matrix_[i][j])) + 1);
            if (main_matrix_[i][j] < 0) n_intd++;
            max_n_intd = (n_intd > max_n_intd) ? n_intd : max_n_intd;
        }
    }
    display_width_ = max_n_intd + pcs + ((pcs > 0) ? 1 : 0);
}

template <size_t D, typename T>
tensor<D, T>& tensor<D, T>::operator= (std::vector< std::vector<T> >& m) {
    int r = m.size();
    int c = 1;
    // iterate through outer list to find largest inner vector
    for (std::vector< std::vector<T> >& x : m) {
        if (x.size() > c) c = x.size();
    }

    if (!Initialisation(r, c)) return *this;
    for (long i = 0; i < r; i++) {
        for (long j = 0; j < m[i].size(); j++) {
            main_matrix_[i][j] = m[i][j];
        }
    }

    return *this;
}

template <size_t D, typename T>
tensor<D, T> tensor<D, T>::operator= (std::initializer_list< std::initializer_list<T> >& il) {
    int r = il.size();
    int c = 1;
    // iterate through outer list to find largest inner list
    for (std::initializer_list< std::initializer_list<T> >& x : il) {
        if (x.size() > c) c = x.size();
    }

    if (!Initialisation(r, c)) return *this;
    typename std::initializer_list< std::initializer_list<T> >::iterator rit = il.begin();
    for (long i = 0; i < r; i++) {
        typename std::initializer_list<T>::iterator cit = rit->begin();
        for (long j = 0; j < c; j++) {
            if (j >= rit->size()) main_matrix_[i][j] = 0;
            else main_matrix_[i][j] = *cit;
            cit++;
        }
        rit++;
    }

    return *this;
}

template <size_t D, typename T>
std::vector<T>& tensor<D, T>::operator[] (const int& index) {
    return main_matrix_[index];
}

template <size_t D, typename T>
void tensor<D, T>::AssignElement(const int& r, const int& c, const T& val) {
    main_matrix_[r][c] = val;
    FlagReset();
}

template <size_t D, typename T>
void tensor<D, T>::Pivot(int k) {
    std::vector< std::vector<double> >& mlu = matrix_lu_;
    std::vector<int>& pr = matrix_pr_;
    std::vector<int>& pc = matrix_pc_;
    T max_val = -128;
    int max_i = k;
    int max_j = k;
    int pk = k;
    for (int i = k; i < n_rows_; i++) {
        for (int j = k; j < n_cols_; j++) {
            if (abs(mlu[pr[i]][pc[j]]) > max_val) {
                max_val = abs(mlu[pr[i]][pc[j]]);
                max_i = i;
                max_j = j;
            }
        }
    }

    pk = pr[k];
    pr[k] = pr[max_i];
    pr[max_i] = pk;
    pk = pc[k];
    pc[k] = pc[max_j];
    pc[max_j] = pk;
}

template <size_t D, typename T>
bool tensor<D, T>::LUDecomposition() {
    std::vector< std::vector<double> >& mlu = this->matrix_lu_;
    std::vector<int>& pr = matrix_pr_;
    std::vector<int>& pc = matrix_pc_;
    int s = std::min(n_rows_, n_cols_);
    double divisor = 1.0;

    if (is_decomposed_ == true) return true;

    for (int i = 0; i < s; i++) {
        if (i <= 0) mlu = main_matrix_;
        else {
            for (int r = i; r < n_rows_; r++) {
                for (int c = i; c < n_cols_; c++) {
                    mlu[pr[r]][pc[c]] -= mlu[pr[r]][pc[i-1]] * mlu[pr[i-1]][pc[c]];
                    if (std::isfinite(mlu[pr[r]][pc[c]]) && (abs(mlu[pr[r]][pc[c]]) < PRECISION)) mlu[pr[r]][pc[c]] = 0;
                    else if (!std::isfinite(mlu[pr[r]][pc[c]])) is_invertible_ = false;
                }
            }
        }

        Pivot(i);
        divisor = 1 / mlu[pr[i]][pc[i]];
        for (int r = i + 1; r < n_rows_; r++) {
            mlu[pr[r]][pc[i]] *= divisor;
            if (std::isfinite(mlu[pr[r]][pc[i]]) && (abs(mlu[pr[r]][pc[i]]) < PRECISION)) mlu[pr[r]][pc[i]] = 0;
            else if (!std::isfinite(mlu[pr[r]][pc[i]])) is_invertible_ = false;
        }
    }

    is_decomposed_ = true;
    return true;
}

template <size_t D, typename T>
static tensor<D, T> tensor<D, T>::LUComposition(const std::vector< std::vector<double> >& m_lu) {
    if (n_rows_ != n_cols_) return NULL;
    tensor<T> ans(n_rows_);
    for (int i = 0; i < n_rows_; i++) {
        for (int j = 0; j < n_cols_; j++) {
            ans.matrix[i][j] = 0;
            for (int k = 0; k <= std::min(i, j); k++) {
                ans.matrix[i][j] += matrix_lu_[(1 + i)*i / 2 + k] * matrix_lu_[(1 + j)*j / 2 + k];
            }
        }
    }
    return ans;
}

template <size_t D, typename T>
void tensor<D, T>::DisplayLU(int pcs) {
    std::vector< std::vector<double> >& mlu = matrix_lu_;
    std::vector<int>& pr = matrix_pr_;
    std::vector<int>& pc = matrix_pc_;
    std::cout.precision(pcs);
    std::cout.setf(std::ios::right);
    std::cout.setf(std::ios::fixed);
    for (int i = 0; i < n_rows_; i++) {
        for (int j = 0; j < n_cols_; j++) {
            std::cout << mlu[pr[i]][pc[j]];
            if (j < n_cols_ - 1) std::cout << '\t';
        }
        std::cout << std::endl;
    }
}

template <size_t D, typename T>
tensor<T> tensor<D, T>::operator*(const tensor &m) {
    if (n_cols_ != m.n_rows_) return NULL;
    tensor<T> ans(n_rows_, m.n_cols_);
    for (int i = 0; i < n_rows_; i++) {
        for (int j = 0; j < m.n_cols_; j++) {
            ans.main_matrix_[i][j] = 0;
            for (int k = 0; k < n_cols_; k++) {
                ans.main_matrix_[i][j] += main_matrix_[i][k] * m.main_matrix_[k][j];
            }
        }
    }
    return ans;
}

template <size_t D, typename T>
double tensor<D, T>::Determinant() {
    std::vector< std::vector<double> >& mlu = matrix_lu_;
    std::vector<int>& pr = matrix_pr_;
    std::vector<int>& pc = matrix_pc_;
    double det = 1;
    if (n_rows_ != n_cols_) return 0;
    for (int i = 0; i < n_rows_; i++) det *= mlu[pr[i]][pc[i]];
    return det;
}

template <size_t D, typename T>
template <size_t dims>
multidimensional_vector<dims, T>::type& tensor<D, T>::select_dimension(multidimensional_vector<dims, T>::type &mv, size_t index) {
    if(dims == 2) return mv[index]
}

long long xinda_math::pow(size_t x, size_t p) {
    if (p == 0) return 1;
    if (p == 1) return x;
    return x * pow(x, p - 1);
}
 */
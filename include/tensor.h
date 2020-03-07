//
// Created by Christopher on 3/4/20.
//

#ifndef TENSOR_TOOLKIT_TENSOR_H
#define TENSOR_TOOLKIT_TENSOR_H


#include <cstddef>
#include <vector>
#include <ostream>
#include <stdexcept>

namespace tensor_toolkit {
    constexpr double PRECISION = 10e-12;
    constexpr unsigned long GIGABYTE = 1073741824;
    constexpr unsigned short NULLPTR_EXC = 0x10;
    constexpr char NULLPTR_EXC_STR[74] = "RUNTIME ERROR: Invalid reference to a tensor element!\r\nProgramme stopped!";
    constexpr unsigned short FAILED_INIT_EXC = 0x10;
    constexpr char FAILED_INIT_STR[64] = "RUNTIME ERROR: Tensor initializatio failed!\r\nProgramme stopped!";
    
    template<typename T>
    class tensor {
    public:
        explicit tensor(std::vector<std::size_t> dimensions);
        
        tensor(std::vector<std::size_t> dimensions, std::vector<T> vec);
        
        class tensor_view {
        public:
            // Constructor
            tensor_view(tensor<T> &vec, std::size_t index, std::size_t dimension);
            
            tensor_view &operator[](std::size_t n_index);
            
            explicit operator T() const;
            
            tensor_view &operator=(T val);
            
            virtual ~tensor_view() = default;
            
            friend std::ostream &operator<<(std::ostream &os, const tensor<T>::tensor_view &tv) {
                if (tv.mValue != nullptr)
                    os << *tv.mValue;
                else {
                    std::cout << NULLPTR_EXC_STR << std::endl;
                    std::exit(NULLPTR_EXC);
                }
                return os;
            };
        
        private:
            tensor<T> &mTensor;
            T *mValue;
            std::size_t mIndex;
            std::size_t mDimension;
        };
        
        tensor<T>::tensor_view operator[](std::size_t index);
        
        void assign_each(const std::vector<T> &vec);
        
        std::vector<T> get_values() const { return mValues; };
    
    private:
        std::vector<T> mValues;
        std::vector<std::size_t> mDimensions;
    };
};

#endif //TENSOR_TOOLKIT_TENSOR_H

// TODO: Need to be defined methods
// tensor<D, T> &operator=(std::vector<std::vector<T> > &m);
// tensor<D, T> &operator=(std::initializer_list<std::initializer_list<T> > &il);
// std::vector<T> &operator[](const int &index);
// void set_element(const int &r = 0, const int &c = 0, const T &val = 0);
// tensor<D, T> operator*(const tensor &m);
// typename multidimensional_vector<D, T>::type &
// select_dimension(typename multidimensional_vector<mD, mT>::type &mv);

//// Overload operators <<
//friend std::ostream &operator<<(std::ostream &stream, tensor<D, T> &vtensor) {
//    stream.setf(std::ios::right);
//    stream.setf(std::ios::fixed);
//    vtensor.SetDisplayWidth(stream.precision());
//    for (int i = 0; i < vtensor.n_rows_; i++) {
//        for (int j = 0; j < vtensor.n_cols_; j++) {
//            stream.width(vtensor.display_width_);
//            stream << vtensor[i][j];
//            if (j < vtensor.n_cols_ - 1) stream << '\t';
//        }
//        stream << std::endl;
//    }
//    return stream;
//};
//
//// Overload operators >>
//friend std::ostream &operator>>(std::ostream &stream, tensor<D, T> &vtensor) {
//    stream.setf(std::ios::right);
//    stream.setf(std::ios::fixed);
//    vtensor.SetDisplayWidth(stream.precision());
//    for (int i = 0; i < vtensor.n_rows_; i++) {
//        for (int j = 0; j < vtensor.n_cols_; j++) {
//            stream.width(vtensor.display_width_);
//            stream << vtensor[i][j];
//            if (j < vtensor.n_cols_ - 1) stream << '\t';
//        }
//        stream << std::endl;
//    }
//    return stream;
//};


// TODO: Realize the following methods when define the child class 'matrix'
// bool LUDecomposition();
// static tensor<D, T> &LUComposition(const std::vector<std::vector<double> > &m_lu = matrix_lu_);
// void DisplayLU(int pcs = 10);
// double Determinant();
// void Pivot(int k);
// std::vector<std::vector<T> > main_matrix_;
// std::vector<std::vector<double> > matrix_lu_;
// std::vector<int> matrix_pr_, matrix_pc_;
// int n_rows_, n_cols_;
// int display_width_;
// bool is_decomposed_;
// bool is_invertible_;
// void FlagReset();
// void SetDisplayWidth(const int &pcs = 2);
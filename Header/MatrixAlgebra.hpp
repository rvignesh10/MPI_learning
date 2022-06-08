#ifndef MATRIX_ALG_HPP
#define MATRIX_ALG_HPP

#include <iostream>
#include "constants.hpp"
/* -------------------------------------------------------------------------------------------------------------------- */
template<typename T, const int N>
class Vector{
    private:
        T m_Vector[N][1];
    public:
        Vector(); // Constructor
        int getLength_(){
            return N;
        }
        void setValue(int idx, T value);
        T getValue(int idx);
        T dotProduct(Vector<T,N> a);
        void displayVector();

};

template<typename T, const int N>
Vector<T,N>::Vector(){
    for (int i=0; i<N; i++){
        setValue(i,0);
    }
    //std::cout << "Vector Created with all Zero values" << std::endl;
}

template<typename T, const int N>
void Vector<T,N>::setValue(int idx, T value){
    if(idx >= N){ std::cerr << "Index length greater than N" << std::endl;}
    else{ m_Vector[idx][0] = value; }
}

template<typename T, const int N> 
T Vector<T,N>::getValue(int idx){
    if(idx>=N){std::cerr << "Index youre trying to access exceeds vector length" << std::endl; return (T)0;}
    else{ return m_Vector[idx][0]; }
}

template<typename T, const int N>
T Vector<T,N>::dotProduct(Vector<T,N> a){
    T sum= (T)0;
    for (int i=0; i<N; i++){
        sum += a.getValue(i)*getValue(i);
    }
    return sum;
}

template<typename T, const int N>
void Vector<T,N>::displayVector(){
    for(int i=0; i<N; i++){std::cout<<m_Vector[i][0]<<std::endl;}
}

/* ---------------------------------------------------------------------------------------------------------------------  */
template<typename T, int m, int n>
class Matrix{
    private:
        T m_Matrix[m][n];
    public:
        Matrix(); // Constructor
        void getSize_(){
            std::cout << "The size is:" << m << "and" << n << std::endl;
        }
};

template<typename T, int m, int n>
Matrix<T,m,n>::Matrix(){
    std::cout << "Matrix created with 0 entries" << std::endl;
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            m_Matrix[i][j] = (T) 0;
        }
    }
}



#endif
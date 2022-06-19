#ifndef MATRIX_ALG_HPP
#define MATRIX_ALG_HPP

#include <iostream>
#include <vector>
#include "constants.hpp"
/* -------------------------------------------------------------------------------------------------------------------- */
template<typename T, const int N>
class Vector{
    private:
        T m_Vector[N];
    public:
        Vector(); // Constructor
        Vector(T *v, int Size); // Constructor to initialize with vales
        void init_();
        int getLength_(){
            return N;
        }
        void setValue(int idx, T value);
        void setVector(T *vec, int Size);
        T getValue(int idx);
        void MultScalar(T scale);
        T dotProduct(Vector<T,N> a);
        T *getData();
        void displayVector();

};

template<typename T, const int N>
Vector<T,N>::Vector(){
   init_();
}

template<typename T, const int N>
Vector<T,N>::Vector(T *v, int Size){
    if(N > Size){
        init_();
        std::cerr << "All values wont be initialized since size of vector > selected vector size" << std::endl;
    }
    else{
        for(int i=0; i<N; i++){
            m_Vector[i] = v[i];
        }
    }
}

template<typename T, const int N>
void Vector<T,N>::init_(){
    for (int i=0; i<N; i++){
        setValue(i,0);
    }
}

template<typename T, const int N>
void Vector<T,N>::setValue(int idx, T value){
    if(idx >= N){ std::cerr << "Index length greater than N" << std::endl;}
    else{ m_Vector[idx] = value; }
}

template<typename T, const int N>
void Vector<T,N>::setVector(T *vec, int Size){
    if(N > Size){
        init_();
        std::cerr << "All values wont be initialized since size of vector > selected vector size" << std::endl;
    }
    else{
        for(int i=0; i<N; i++){
            m_Vector[i] = vec[i];
        }
    }
}

template<typename T, const int N> 
T Vector<T,N>::getValue(int idx){
    if(idx>=N){std::cerr << "Index youre trying to access exceeds vector length" << std::endl; return (T)0.0;}
    else{ return m_Vector[idx]; }
}

template<typename T, const int N>
void Vector<T,N>::MultScalar(T scale){
    for(int i=0; i<N; i++){
        m_Vector[i] *= scale;
    }
}

template<typename T, const int N>
T Vector<T,N>::dotProduct(Vector<T,N> a){
    T sum= (T)0.0;
    for (int i=0; i<N; i++){
        sum += a.getValue(i)*getValue(i);
    }
    return sum;
}


template<typename T, const int N>
T *Vector<T,N>::getData(){
    T *temp_vector = NULL;
    temp_vector = m_Vector;
    return temp_vector;
}

template<typename T, const int N>
void Vector<T,N>::displayVector(){
    for(int i=0; i<N; i++){std::cout<<m_Vector[i]<<" ";}
}

/* ---------------------------------------------------------------------------------------------------------------------  */
template<typename T, const int m, const int n>
class Matrix{
    private:
        T m_Matrix[m][n];
        bool SQUARE;
    public:
        Matrix(); // Constructor
        Matrix(T *rowMajor, int Size); // imports matrix data from a rowMajor format
        void getSize_(){
            std::cout << "The size is:" << m << " x " << n << std::endl;
        }
        Matrix(int i); // Matrix to generate i*Identity matrix 
        void init_();
        void setValue(int i, int j, T val);
        void setRow(int row_idx, Vector<T,n> r);
        void setColumn(int col_idx, Vector<T,m> c);
        void getRow(int row_idx, Vector<T,n> &r);
        void getColumn(int col_idx, Vector<T,m> &c);
        void Transpose(Matrix<T,n,m> &m_Matrix_T);
        void Add(Matrix<T,m,n> M_in, Matrix<T,m,n> &M_out);
        void Subtract(Matrix<T,m,n> M_in, Matrix<T,m,n> &M_out);
        void MultScalar(T scale, Matrix<T,m,n> &M_out);
        void MultVector(Vector<T,n> u, Vector<T,m> &u_out);
        template<const int C>
        void MultMatrix(Matrix<T,n,C> M_in, Matrix<T,m,C> &M_out);
        bool checkSquare(){return SQUARE;}
        T *exportRowMajor();
        void setRowMajor(T *rowMajor, int Size);
};

template<typename T, const int m, const int n>
Matrix<T,m,n>::Matrix(){
    init_();
    
}

template<typename T, const int m, const int n>
Matrix<T,m,n>::Matrix(T *rowMajor, int Size){
    if (m*n != Size){
        std::cerr<<"Matrix size is incompatible" << std::endl;
    }
    else {
        int idx = 0;
        for (int i=0; i<m; i++){
            for (int j=0; j<n; j++){
                m_Matrix[i][j] = rowMajor[idx];
                ++idx;
                }
        }
        SQUARE = (m==n)?true:false;
    }
}

template<typename T, const int m, const int n>
Matrix<T,m,n>::Matrix(int i){
    init_();
    if(SQUARE){
        for(int idx=0; idx<n; idx++){
            m_Matrix[idx][idx] = (T)(i)*(T)(1.0);
        }
    }
    else{
        std::cerr<<"Not a square matrix and hence it is initialized with zeros" << std::endl;
    }
}

template<typename T, const int m, const int n>
void Matrix<T,m,n>::init_(){
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            m_Matrix[i][j] = (T) 0.0;
        }
    }
    SQUARE = (m==n)?true:false;
}

template<typename T, const int m, const int n>
void Matrix<T,m,n>::setValue(int i, int j, T val){
    if(i>=m || j>=n){
        std::cerr << "Index exceeds matrix dimensions" << std::endl;
    }
    else {
        m_Matrix[i][j] = val;
    }
}

template<typename T, const int m, const int n>
void Matrix<T,m,n>::setRow(int row_idx, Vector<T,n> r){
    if (row_idx >= m){
        std::cerr << "Row index exceeds Matrix dimensions" << std::endl;
    }
    else {
        for(int i = 0; i < n; i++){
            m_Matrix[row_idx][i] = r.getValue(i);
        }
    }
}

template<typename T, const int m, const int n>
void Matrix<T,m,n>::setColumn(int col_idx, Vector<T,m> c){
    if (col_idx >= n){
        std::cerr << "Column index exceeds Matrix dimensions" << std::endl;
    }
    else {
        for (int i = 0; i < m; i++){
            m_Matrix[i][col_idx] = c.getValue(i);
        }
    }
}

template<typename T, const int m, const int n>
void Matrix<T,m,n>::getRow(int row_idx, Vector<T,n> &r){
    if(row_idx >= m){
        std::cerr << "Row index exceeds matrix dimensions" << std::endl;
    }
    else {
        for(int i=0; i<n; i++){
            r.setValue(i,m_Matrix[row_idx][i]);
        }
    }
}

template<typename T, const int m, const int n>
void Matrix<T,m,n>::getColumn(int col_idx, Vector<T,m> &c){
    if(col_idx >= n){
        std::cerr << "Column index exceeds matrix dimensions" << std::endl;
    }
    else {
        for(int i=0; i<m; i++){
            c.setValue(i,m_Matrix[i][col_idx]);
        }
    }
}

template<typename T, const int m, const int n>
void Matrix<T,m,n>::Transpose(Matrix<T,n,m> &m_Matrix_T){
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            m_Matrix_T.setValue(j,i,m_Matrix[i][j]);
        }
    }
}

template<typename T, const int m, const int n>
void Matrix<T,m,n>::Add(Matrix<T,m,n> M_in, Matrix<T,m,n> &M_out){
    T *rm1 = exportRowMajor();
    T *rm2 = M_in.exportRowMajor();
    T *vec = (T*)malloc((m*n)*sizeof(T));
    for (int i=0; i < m*n; i++){
        vec[i] = rm1[i] + rm2[i];
    }
    Matrix<T,m,n> test(vec, m*n);
    //M_out.setRowMajor(vec, m*n);
    M_out = test;
}

template<typename T, const int m, const int n>
void Matrix<T,m,n>::Subtract(Matrix<T,m,n> M_in, Matrix<T,m,n> &M_out){
    T *rm1 = exportRowMajor();
    T *rm2 = M_in.exportRowMajor();
    T *vec = (T*)malloc((m*n)*sizeof(T));
    for (int i=0; i < m*n; i++){
        vec[i] = rm1[i] - rm2[i];
    }
    Matrix<T,m,n> test(vec);
    //M_out.setRowMajor(vec, m*n);
    M_out = test;
}

template<typename T, const int m, const int n>
void Matrix<T,m,n>::MultScalar(T scale, Matrix<T,m,n> &M_out){
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            M_out.setValue(i,j,m_Matrix[i][j]*scale);
        }
    }
}

template<typename T, const int m, const int n>
void Matrix<T,m,n>::MultVector(Vector<T,n> u, Vector<T,m> &u_out){
    for(int i=0; i<m; i++){
        Vector<T,n> row ;
        getRow(i,row);
        T val = u.dotProduct(row);
        u_out.setValue(i,val);
    }
}

template<typename T, const int m, const int n>
template<const int C>
void Matrix<T,m,n>::MultMatrix(Matrix<T,n,C> M_in, Matrix<T,m,C> &M_out){
    for(int i=0; i<m; i++){
        Vector<T,n> row;
        getRow(i,row);
        for(int j=0; j<C; j++){
            Vector<T,n> column;
            M_in.getColumn(j,column);
            T val = row.dotProduct(column);
            M_out.setValue(i,j,val);
        }
    }
}

/// exports the matrix as a pointer to an array. The rows of the matrix are organized one next to each other.
template<typename T, const int m, const int n>
T *Matrix<T,m,n>::exportRowMajor(){
    T *rowMajor = (T*)malloc((m*n)*sizeof(T));
    int idx = 0;
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            rowMajor[idx] = m_Matrix[i][j];
            ++idx;
        }
    }
    return rowMajor;
}

template<typename T, const int m, const int n>
void Matrix<T,m,n>::setRowMajor(T *rowMajor, int Size){
    if(m*n != Size){
        std::cerr << "Array size is not compatible with the matrix for allocation" << std::endl;
    }
    else {
        int idx = 0;
        for(int i=0; i<m; i++){
            for(int j=0; j<n; j++){
                m_Matrix[i][j] = rowMajor[idx];
                ++idx;
            }
        }
    }
}

#endif
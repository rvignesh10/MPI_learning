#include <iostream>
#include <string>
#include "../Header/MatrixAlgebra.hpp"

template<typename T, const int m, const int n>
void printMat(Matrix<T,m,n> M, std::string S){
    std::cout << "This Matrix after, " << S << " operation looks like: " << std::endl; 
    for (int i=0; i<m; i++){
        Vector<T,n> r;
        M.getRow(i,r);
        r.displayVector();
        std::cout << std::endl;
    }    
}

int main(){
    const int N = 9;
    const int n = 3;
    double row[N];
    for(int i=0; i<N; i++){
        row[i] = 1.0*((double)(i)+1.0);
    }
    // Test initialization
    Matrix<double,n,n> m(row,N);
    printMat(m,"print");

    // Test Transpose
    Matrix<double,n,n> testMat;
    m.Transpose(testMat);
    printMat(testMat,"Transpose");

    // Test Addition and Subtraction
    testMat.Add(m,testMat);
    printMat(testMat, "Addition");

    // Test Vector Multiplication
    Vector<double,n> u(row,N);
    m.MultVector(u,u);
    std::cout<< "The matrix multiplication of a vector operation is: " << std::endl;
    u.displayVector();
    std::cout<<std::endl;

    // Test identity matrix
    Matrix<double,n,n> eye(1);
    printMat(eye,"Identity");

    // Test Matrix Multiplication - scalar
    eye.MultScalar(4.0,eye);
    printMat(eye,"Matrix scalar multiplication");

    // Test Matrix vector multiplication
    u.init_();
    u.setVector(row,N);
    eye.MultVector(u,u);
    u.displayVector();
    std::cout<<std::endl;

    // Test Matrix multiplication 
    testMat.MultMatrix(m,testMat);
    printMat(testMat,"Matrix multiplication");

    // Test exporting as rowMajor
    return 0;
}
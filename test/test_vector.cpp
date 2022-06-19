#include <iostream>
#include <vector>
#include "../Header/MatrixAlgebra.hpp"

int main(){
    const int N = 4;
    Vector<double,N> v;

    double *stdVector = v.getData();

    for(int i=0; i<N; i++){
        std::cout<<stdVector[i] << " " ;
    }
    std::cout<<std::endl;

}
//#include "../Header/MatrixAlgebra.hpp"
#include "../Header/IntegrationRule.hpp"
#include <iostream>

int main(){
    Vector<float,INT_PTS+TWO> vec;
    Vector<float,INT_PTS+TWO> wts;
    IntegrationRule<float,INT_PTS+TWO>(-0.1f,0.1f,vec,wts);
    vec.displayVector();
    std::cout<< std::endl;
    wts.displayVector();
    return 0;
}
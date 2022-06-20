#include <iostream>
#include "../Header/MatrixAlgebra.hpp"
#include "../Header/FiniteElement.hpp"

int main(){
    const int N = 5;
    Mesh<N> m(-10.0,10.0);
    m.gen_1D_mesh();
    
    const int order = 2;
    FiniteElementSpace<order,N> fes(m);
    for (int i=0; i<N; i++){
        fes.display_FE_Info(i);
        std::cout << std::endl;
    }
    
    return 0;
}
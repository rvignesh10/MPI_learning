#include <iostream>
#include <mpi.h>
#include "../Header/MatrixAlgebra.hpp"
#include "../Header/FiniteElement.hpp"

int main(){

    const int order = 3;
    {   // check serial mode initialization
        Mesh<FIVE*TWO> m(-10.0,10.0);
        FiniteElementSpace<order,FIVE*TWO> fes(m);
        for (int i=0; i<FIVE*TWO; i++){
            fes.display_FE_Info(i);
            std::cout << std::endl;
        }
    }
    // check parallel mode initialization
    // MPI Mesh testing
    MPI_Init(NULL,NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int N = 20; // assume it uses 4 processors

    double g_xlim1 = -10.0;
    double g_xlim2 = 10.0;

    // Finding processors local domain
    double DX = (g_xlim2-g_xlim1)/(double)world_size;
    double l_xlim1 = g_xlim1 + (double)world_rank*DX;
    double l_xlim2 = l_xlim1 + DX;

    // Finding number of elements each domain has
    const int N_local = 5;
    Mesh<N_local> m(l_xlim1,l_xlim2);
    FiniteElementSpace<order,N_local> fes(m,world_rank);

    if (world_rank == 2){
        std::cout << "World's rank is: " << world_rank << std::endl;
        fes.display_FE_Info(0);
    }

    return 0;
}
#include <iostream>
#include <mpi.h>
#include <string>
#include "../../Header/FiniteElement.hpp"
using namespace std;

int main(int argc, char **argv){

    const int numProcesses = 2;
    // if (argc >=3){
    //     numProcesses = stoi(argv[2]);
    // }
    // MPI Mesh testing
    MPI_Init(NULL,NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    double g_xlim1 = -1.0;
    double g_xlim2 = 1.0;

    // Finding processors local domain
    double DX = (g_xlim2-g_xlim1)/(double)world_size;
    double l_xlim1 = g_xlim1 + (double)world_rank*DX;
    double l_xlim2 = l_xlim1 + DX;

    // Finding number of elements each domain has
    const int N_local = 2;
    const int order = 2;

    Mesh<N_local> m(l_xlim1,l_xlim2);
    FiniteElementSpace<order,N_local> fes(m,world_rank);
    //fes.display_FE_Info(1);

    // Assign Global Matrix for Diffusion Integration
    Matrix<double,(numProcesses*N_local+ONE)+(numProcesses*N_local)*(order-ONE),
           (numProcesses*N_local+ONE)+(numProcesses*N_local)*(order-ONE)> DiffusionMat;

    AddDomainIntegrators Integrator;
    AppendList *diff = new AppendList;
    //std::cout << "You are in main where diff is: " << diff << std::endl;
    Integrator.DiffusionIntegrator(fes,&diff);
    //std::cout << "You came out and diff is: " << diff << std::endl;
    DiffusionMat.AssignFromList(diff);
    //DiffusionMat.displayMatrix();
    delete diff;
}
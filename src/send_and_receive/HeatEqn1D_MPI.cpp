#include <iostream>
#include <mpi.h>
#include <math.h>
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

    double nu = 1.0;
    // Domain set up
    double g_xlim1 = -1.0;
    double g_xlim2 = 1.0;
    double tlim1   = 0.0;
    double tfinal  = 0.02;
    double r       = 0.5;

    // Finding processors local domain
    double DX = (g_xlim2-g_xlim1)/(double)world_size;
    double l_xlim1 = g_xlim1 + (double)world_rank*DX;
    double l_xlim2 = l_xlim1 + DX;

    // Finding number of elements each domain has
    const int N_local = 2;
    const int order = 1;
    const int NTot  = (numProcesses*N_local+ONE)+(numProcesses*N_local)*(order-ONE);

    Mesh<N_local> m(l_xlim1,l_xlim2);
    double dx;
    m.getDX(dx);
    double dt = r*pow(dx,2.0)/nu;

    FiniteElementSpace<order,N_local> fes(m,world_rank);

    // Assign Global Matrix for Diffusion Integration
    Matrix<double,NTot,NTot> DiffusionMat;
    Matrix<double,NTot,NTot> MassMat;
    

    AddDomainIntegrators Integrator;
    AppendList *diff ;
    AppendList *mass ;
    Integrator.DiffusionIntegrator(fes,&diff);
    Integrator.MassIntegrator(fes,&mass);

    DiffusionMat.AssignFromList(diff);
    MassMat.AssignFromList(mass);

    double *diff_rm = (double*)malloc((NTot*NTot)*sizeof(double));
    double *mass_rm = (double*)malloc((NTot*NTot)*sizeof(double));

    // define global row major form of matrices
    double *g_diff_rm = (double*)malloc((NTot*NTot)*sizeof(double));
    double *g_mass_rm = (double*)malloc((NTot*NTot)*sizeof(double));

    diff_rm = DiffusionMat.exportRowMajor();
    mass_rm = MassMat.exportRowMajor();
    
    // This piece sends out the matrix as a rowMajor pointer array and adds it
    // with the local diffusion and mass matrices (also stored as rowMajor pointer array).
    // Then it converts it to a global diffusion/mass matrix. 
    if(world_size>1){
        if(world_rank==0){
            for (int j=0; j<NTot*NTot; j++){
                g_diff_rm[j] = diff_rm[j];
                g_mass_rm[j] = mass_rm[j];
            }
            std::cout << std::endl;
            for(int i=1; i<world_size; i++){
                MPI_Recv(diff_rm,NTot*NTot,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Recv(mass_rm,NTot*NTot,MPI_DOUBLE,i,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

            for (int j=0; j<NTot*NTot; j++){
                g_diff_rm[j] += diff_rm[j];
                g_mass_rm[j] += mass_rm[j];
            }
    
            }
            Matrix<double,NTot,NTot> g_diffMat(g_diff_rm,NTot*NTot);
            Matrix<double,NTot,NTot> g_massMat(g_mass_rm,NTot*NTot);

        }
        else {
            MPI_Send(diff_rm,NTot*NTot,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
            MPI_Send(mass_rm,NTot*NTot,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
        }
    }
    
    free(g_diff_rm);
    free(g_mass_rm);
    free(diff_rm);
    free(mass_rm);
    delete diff;
    delete mass;
}
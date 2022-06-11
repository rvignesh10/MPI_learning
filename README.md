# MPI_learning
This repository is to hold files and instructions on how to use MPI for Numerical Solvers. *Enter*
Please make sure you have MPICH-4.0.2 installed and you can find the package at https://www.mpich.org

## Installation
> 1. Download package and unzip it. 
> 2. Perform ``` ./configure --prefix=/installation/directory/path --disable-fortran ```
> 3. Run ``` make; sudo install ```
> 4. Check if you've installed correctly using ``` mpiexec --version ```
> 5. Ensure that you have set the compiler to use ``` export MPICXX=~/bin/mpicxx ```

## Header
This folder will have the required header files that needs to be included in order to run test and *Enter*
example cases. The header files included are for 1D mesh generation, FiniteElementSpace creation and *Enter*
IntegrationRule for these discretized spaces.

## test
The test folder will help in testing some of the changes added in the header files and make sure the *Enter*
the functionality of the defined classes are intact. 

## src
This will contain the source code for the examples we will use to learn about MPI and how to use it for *Enter*
Numerical Methods and processes. If you either make a change to the source code or add a new source code *Enter*
you can compile and run the executable typing the following in your command window prompt *Enter*
``` mpicxx <file_name.cpp> -o <executable_name> ```
``` mpirun -np x ./executable_name ```

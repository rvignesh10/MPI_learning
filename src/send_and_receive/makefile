EXECS=NumericalIntegration HeatEqn1D_MPI
MPICXX?=mpicxx

all: ${EXECS}

NumericalIntegration: NumericalIntegration.cpp
	${MPICXX} -o NumericalIntegration NumericalIntegration.cpp

HeatEqn1D_MPI: HeatEqn1D_MPI.cpp
	${MPICXX} -o HeatEqn1D_MPI HeatEqn1D_MPI.cpp

clean:
	rm -f ${EXECS}

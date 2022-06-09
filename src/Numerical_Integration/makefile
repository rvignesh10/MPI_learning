EXECS=NumericalIntegration
MPICXX?=mpicxx

all: ${EXECS}

NumericalIntegration: NumericalIntegration.cpp
	${MPICXX} -o NumericalIntegration NumericalIntegration.cpp

clean:
	rm -f ${EXECS}

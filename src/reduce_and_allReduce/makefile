EXECS=NumericalIntegration_allReduce
MPICXX?=mpicxx

all: ${EXECS}

NumericalIntegration_allReduce: NumericalIntegration_allReduce.cpp
	${MPICXX} -o NumericalIntegration_allReduce NumericalIntegration_allReduce.cpp

clean:
	rm -f ${EXECS}

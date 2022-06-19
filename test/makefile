EXECS=Integration_test mesh_test test test_vector test_matrix
MPICXX?=mpicxx

all: ${EXECS}

Integration_test: Integration_test.cpp
	${MPICXX} -o Integration_test Integration_test.cpp

mesh_test: mesh_test.cpp
	${MPICXX} -o mesh_test mesh_test.cpp

test: test.cpp
	${MPICXX} -o test test.cpp

test_vector: test_vector.cpp
	${MPICXX} -o test_vector test_vector.cpp

test_matrix: test_matrix.cpp
	${MPICXX} -o test_matrix test_matrix.cpp

clean:
	rm -f ${EXECS}

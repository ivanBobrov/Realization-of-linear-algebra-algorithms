#ifndef JACOBI_H
#define JACOBI_H

#include <mpi.h>
#include "define.h"
#include "vector.h"
#include "sym_matrix.h"
#include "csr_matrix.h"


vector *	 jacobi_solve ( sym_matrix *, vector *, double);
vector *	 jacobi_solve ( csr_matrix *, vector *, double);
vector * mpi_jacobi_solve ( csr_matrix *, vector *, double, int procNumber, int currentProc );
vector *	 jacobi_iteration ( sym_matrix * A, vector * f, vector * prev, vector * next );
vector *	 jacobi_iteration ( csr_matrix * A, vector * f, vector * prev, vector * next );
vector * mpi_jacobi_iteration ( csr_matrix * A, vector * f, vector * prev, vector * next, int procNumber, int currentProc );
double calculateNextComponent ( csr_matrix * A, vector * f, vector * prev, int component );
double get_norm ( vector * next, vector * prev );

#ifdef GCC
double fabs ( double );
#endif

#endif
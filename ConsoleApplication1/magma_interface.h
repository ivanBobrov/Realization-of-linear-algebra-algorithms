#ifndef MAGMA_INTERFACE_H
#define MAGMA_INTERFACE_H

#include "vector.h"
#include "csr_matrix.h"

typedef int magma_int_t;

extern "C" magma_int_t magma_sgesv( magma_int_t n, magma_int_t nrhs
								  , float *A, magma_int_t lda
								  , magma_int_t *ipiv 
								  , float *B, magma_int_t ldb
								  , magma_int_t *info );

vector * magma_jacobi_solve ( csr_matrix *, vector * );
vector * magma_jacobi_solve ( float * M, float * V, int size ); 
float * packMatrixToBLAS ( csr_matrix * input );
float * packVectorToBLAS ( vector * input );
vector * unpackVectorFromBLAS ( float * input, int size );
void deleteBlasMatrix ( float * R );
void deleteBlasVector ( float * input );

#endif
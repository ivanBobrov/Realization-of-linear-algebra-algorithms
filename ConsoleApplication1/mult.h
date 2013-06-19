#ifndef MULT_H
#define MULT_H

#include <cuda.h>
#include "csr_matrix.h"
#include "vector.h"
#include "define.h"

vector * cuda_jacobi_solve ( csr_matrix *, vector *, double);

#endif
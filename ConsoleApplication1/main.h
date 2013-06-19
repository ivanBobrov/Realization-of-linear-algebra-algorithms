#ifndef MAIN_H
#define MAIN_H

#define DEFAULT_DIMENSION 10
#define MAX_TEST_DIMENSION 1025
#define DIMENSION_MULTIPLIER 2

#define MAX_MODULE_NUMBER 100
#define SINGULARITY_NUMBER 100
#define JACOBI_EPSILON 0.000001

#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <Windows.h>

#include "vector.h"
#include "sym_matrix.h"
#include "jacobi.h"
#include "csr_matrix.h"
#include "define.h"
#include "mult.h"
#include "magma_interface.h"

int PROC_NUMBER;
int CURRENT_PROCESS;
int DIMENSION = DEFAULT_DIMENSION;

csr_matrix * A;
vector * result;
vector * check;
vector * b;
double timeStart, timeEnd;

void printResult( void );
csr_matrix * loadMatrixA ( void );
vector * loadVectorB ( void );
void memoryRelease ( void );
void printTimeToFile( char * filename );

double * packMatrix ( csr_matrix * input );
csr_matrix * unPackMatrix ( double * input );

void sendMatrixToAllProcesses( csr_matrix ** );
void sendVectorToAllProcesses( vector ** );

int bufSize ( double * array, csr_matrix * matrix );

#endif
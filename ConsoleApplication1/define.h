#ifndef DEFINE_H
#define DEFINE_H

//#define MPI
#define CUDA
//#define MAGMA
#define TEST

//#define VERBOSE
//#define MPI_VERBOSE
//#define CUDA_VERBOSE
#define CUDA_INFO
//#define CUDA_DEBUG
#define PRINT_RESULT

#define GENERATE_MATRIX
#define GENERATE_VECTOR

//#define GCC

#define ROOT_PROCESS 0
#define IS_ROOT_PROCESS CURRENT_PROCESS == ROOT_PROCESS
#define IS_WORK_PROCESS CURRENT_PROCESS != ROOT_PROCESS
#define IS_SUCCESSIVE PROC_NUMBER == 1

#define THREADS_PER_BLOCK 512

#endif
#include "mult.h"

void printCudaInfo ( void );
double cuda_get_norm ( float * next, float * prev, int size );
float * packMatrixRowA ( csr_matrix * , int rowToPack );
float * packVector ( vector * );
float * packMatrixA ( csr_matrix * );
vector * unpackVector ( float * input, int size );

float * cuda_jacobi_iteration ( float * A, float * f, float * prev, int size );
 __global__ void jacobiNextComponent ( float * matrix, float * f, float * prev, float * next, int size );


//==============================================================
//============== CUDA операции якоби ===========================
//==============================================================

//============== Kernel ========================================
__global__ void jacobiNextComponent ( float * matrix, float * f, float * prev, float * next, int size ) {
	int i_block = blockIdx.x * blockDim.x + threadIdx.x;
	int component = blockIdx.y;
	float matrixElement;
	float vectorElement;
	__shared__ float R;
	__shared__ float multiplicatedVector[THREADS_PER_BLOCK];
	//__shared__ float device_next[THREADS_PER_BLOCK];

	if ( i_block < size ) {
		matrixElement = matrix[component*size + i_block];
		vectorElement = prev[i_block];

		multiplicatedVector[i_block] = matrixElement * vectorElement;
	}
	if ( i_block == 0 ) {
		R = -f[component];
		
		for ( int j = 0; j < size; j++ )
			R += multiplicatedVector[j];
		
		R -= matrix[component*size + component] * prev[component];
		R /= - matrix[component*size + component];
		
		next[component] = R;
	}
}

//==============================================================
//============= System solution ================================
vector * cuda_jacobi_solve ( float * A, float * f, int size, double eps) {
	float * prev = new float[size];
	float * next = NULL;
	for ( int i = 0; i < size; i++ )
		prev[i] = f[i];			// «атравка
	
	double norm = eps;

#ifdef CUDA_INFO
	printCudaInfo();
#endif

	while ( norm >= eps ) {

		if (next != NULL) {
			delete next;
			next = NULL;
		}
		
		next = cuda_jacobi_iteration( A, f, prev, size );
		norm = cuda_get_norm (next, prev, size);
		for ( int i = 0; i < size; i++ )
			prev[i] = next[i];
		
		#ifdef CUDA_VERBOSE
		printf( "norm = %f\n", norm );
		#endif
	}
	
	vector * result = unpackVector(prev, size);
	if (prev != NULL) delete[] prev;
	prev = NULL;
	if (next != NULL) delete[] next;
	next = NULL;

	return result;
}

//==============================================================
//============ Single iteration ================================
float * cuda_jacobi_iteration ( float * A, float * f, float * prev, int size ) {
#ifdef CUDA_VERBOSE
	printf("CUDA iteration started\n");
#endif
	
	dim3 threads;
	dim3 blocks;
	
	if ( size < THREADS_PER_BLOCK )
		 threads = dim3 ( size );
	else 
		threads = dim3 ( THREADS_PER_BLOCK );
	
	blocks  = dim3 ( (unsigned int)ceil( ((double)size) / ((double)THREADS_PER_BLOCK) ), size );

#ifdef CUDA_VERBOSE
	printf("registered threads : [ %d x %d ]\nregistered blocks  : [ %d x %d ]\n", threads.x, threads.y, blocks.x, blocks.y );
#endif
	
	float * hostMatrixA = A;
	float * hostVectorF = f;
	float * hostVectorPrev = prev;
	float * hostVectorNext = new float[size];
	
	float * deviceMatrixA;
	float * deviceVectorF;
	float * deviceVectorPrev;
	float * deviceVectorNext;
	
	cudaMalloc( (void **)&deviceMatrixA, sizeof(float)*size*size );
	cudaMalloc( (void **)&deviceVectorF, sizeof(float)*size );
	cudaMalloc( (void **)&deviceVectorPrev, sizeof(float)*size );
	cudaMalloc( (void **)&deviceVectorNext, sizeof(float)*size );
	
	cudaMemcpy( deviceMatrixA
			  , hostMatrixA
			  , sizeof(float)*size*size
			  , cudaMemcpyHostToDevice );

	cudaMemcpy( deviceVectorF
			  , hostVectorF
			  , sizeof(float)*size
			  , cudaMemcpyHostToDevice );

	cudaMemcpy( deviceVectorPrev
			  , hostVectorPrev
			  , sizeof(float)*size
			  , cudaMemcpyHostToDevice );
	
	jacobiNextComponent<<<blocks, threads>>>( deviceMatrixA
											, deviceVectorF
											, deviceVectorPrev
											, deviceVectorNext
											, size );
	
	cudaEvent_t syncEvent;
	cudaEventCreate(&syncEvent);
	cudaEventRecord(syncEvent, 0);
	cudaEventSynchronize(syncEvent);
	
	cudaMemcpy( hostVectorNext
			  , deviceVectorNext
			  , sizeof(float)*size
			  , cudaMemcpyDeviceToHost );
	
#ifdef CUDA_DEBUG
	for ( int i = 0; i < size; i++ )
		printf("next[%d]: %f\n", i, hostVectorNext[i]);
#endif
	
	cudaFree(deviceMatrixA);
	cudaFree(deviceVectorF);
	cudaFree(deviceVectorPrev);
	cudaFree(deviceVectorNext);
	
#ifdef CUDA_VERBOSE
	printf("CUDA iteration ended\n");
#endif

	return hostVectorNext;
}

//==============================================================
//================ Common functions ============================
//==============================================================

//================ Norm of residual ============================
double cuda_get_norm ( float * next, float * prev, int size ) {
	// —читаем норму, как максимальную поэлементную разность
	double result;
	result = fabs ( prev[0] - next[0] );
	
	for ( int i = 0; i < size; i++ ) {
		if ( fabs ( prev[i] - next[i] ) > result )
			result = fabs ( prev[i] - next[i] );
	}

	return result;
}

//==============================================================
//================ Packing objects to arrays ===================

float * packMatrixRowA ( csr_matrix * A, int rowToPack ) {
	
	float * result = new float[A->getColCount()];
	for ( int i = 0; i < A->getColCount(); i++ )
		result[i] = (float)A->getElement(rowToPack, i);

	return result;
}

//==============================================================

float * packVector ( vector * vect ) {
	
	float * result = new float[vect->N];
	for ( int i = 0; i < vect->N; i++ )
		result[i] = (float)vect->P[i];

	return result;
}

//==============================================================

float * packMatrixA ( csr_matrix * input ) {
	int size = input->getColCount()*input->getRowCount();
	float * result = new float[size];
	
	for ( int i = 0; i < input->getRowCount(); i++ )
		for ( int j = 0; j < input->getColCount(); j++ )
			result[i*input->getColCount() + j] = (float)(input->getElement(i,j));

	return result;
}

//==============================================================

vector * unpackVector ( float * input, int size ) {
	vector * result = new vector(size);
	
	for ( int i = 0; i < size; i++ )
		result->P[i] = input[i];
	
	return result;
}

//==============================================================
//============= Print cuda support information =================
void printCudaInfo ( void ) {
	int deviceCount;
	cudaError_t err;
	err = cudaGetDeviceCount ( &deviceCount );
	if ( deviceCount <= 0 ) {
		printf("ERROR: No CUDA device found! ( deviceCount = %d ) Error : %s\n", deviceCount, cudaGetErrorString(err));
		return;
	} else
		printf("Found %d devices\n", deviceCount);

	for ( int i = 0; i < deviceCount; i++ ) {
		cudaDeviceProp devProp;
		cudaGetDeviceProperties ( &devProp, i );
	
		printf("CUDA info :\n");
		printf("   Compute capability      : %d.%d\n", devProp.major, devProp.minor );
		printf("   Device name             : %s\n", devProp.name );
		printf("   Total global memory     : %d\n", devProp.totalGlobalMem );
		printf("   Shared memory per block : %d\n", devProp.sharedMemPerBlock );
		printf("   Max threads per block   : %d\n", devProp.maxThreadsPerBlock );
		printf("   Total constant memory   : %d\n", devProp.totalConstMem );
		printf("   MaxGridSize             : [ %d x %d x %d ]\n", devProp.maxGridSize[0], devProp.maxGridSize[1], devProp.maxGridSize[2] );
		printf("=================================\n");
	}
}

//==============================================================

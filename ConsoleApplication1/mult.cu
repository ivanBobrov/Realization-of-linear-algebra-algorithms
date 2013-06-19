#include "mult.h"

void printCudaInfo ( void );
double cuda_get_norm ( vector * next, vector * prev );
float * packMatrixRowA ( csr_matrix * , int rowToPack );
float * packVector ( vector * );
float * packMatrixA ( csr_matrix * );
vector * unpackVector ( float * input, int size );

vector * cuda_jacobi_iteration ( csr_matrix * A, vector * f, vector * prev );
double calculateNextComponentCuda ( csr_matrix * A, vector * f, vector * prev, int component );
 __global__ void matrVectMult( float * row, float * column, float * result, int size );
 __global__ void jacobiNextComponent ( float * matrix, float * f, float * prev, float * next, int size );

//==========================================================
//================= MPI CUDA Jacobi ========================
//==========================================================

// ================ Kernel =================================
 __global__ void matrVectMult( float * row, float * column, float * result, int size ) {

	 if ( (threadIdx.x < size) && (threadIdx.y == 0) )
		result[threadIdx.x] = row[threadIdx.x] * column[threadIdx.x];
	
}

// ================ Kernel launch ==========================

double calculateNextComponentCuda ( csr_matrix * A, vector * f, vector * prev, int component ) {
	int size = f->N;

	float * hostNextComponent = new float[size];
	float * hostPackedMatrixRowA = packMatrixRowA( A, component );
	float * hostPackedVectorPreviousValue = packVector(prev);

	float * deviceNextComponent = NULL;
	float * devicePackedMatrixRowA;
	float * devicePackedVectorPreviousValue;

	cudaMalloc( (void **)&devicePackedMatrixRowA, sizeof(float)*size );
	cudaMalloc( (void **)&devicePackedVectorPreviousValue, sizeof(float)*size );
	cudaMalloc( (void **)&deviceNextComponent, sizeof(float)*size );
	 
	cudaMemcpy( devicePackedMatrixRowA
			  , hostPackedMatrixRowA
			  , sizeof(float)*size
			  , cudaMemcpyHostToDevice );
	
	cudaMemcpy( devicePackedVectorPreviousValue
			  , hostPackedVectorPreviousValue
			  , sizeof(float)*size
			  , cudaMemcpyHostToDevice );

	dim3 grid((size+255)/256, 1, 1);
	dim3 threads(256, 1, 1);

	matrVectMult<<<1, 100>>>( devicePackedMatrixRowA
									 , devicePackedVectorPreviousValue
									 , deviceNextComponent
									 , size );
	
	cudaEvent_t syncEvent;

	cudaEventCreate(&syncEvent);
	cudaEventRecord(syncEvent, 0);
	cudaEventSynchronize(syncEvent);

	cudaMemcpy( hostNextComponent
			  , deviceNextComponent
			  , sizeof(float)*size
			  , cudaMemcpyDeviceToHost );

	
	double result = (double)hostNextComponent[0];
	
#ifdef CUDA_VERBOSE
	for ( int i = 0; i < size; i++ )
		if ( (hostNextComponent[i] > 100) || (hostNextComponent[i] < 0) )  
			printf("CUDA result: %f -- vs -- %f\n", hostNextComponent[i], hostPackedMatrixRowA[i]);
	printf("=======================================\n");
#endif

	cudaEventDestroy(syncEvent);
	cudaFree(deviceNextComponent);
	cudaFree(devicePackedMatrixRowA);
	cudaFree(devicePackedVectorPreviousValue);
	
	delete hostNextComponent;
	delete hostPackedMatrixRowA;
	delete hostPackedVectorPreviousValue;
	
	return result;
}


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
vector * cuda_jacobi_solve ( csr_matrix * A, vector * f, double eps) {
	vector * prev = new vector(*f);
	vector * next = new vector( f->N );
	double norm = eps;

#ifdef CUDA_INFO
	printCudaInfo();
#endif

	while ( norm >= eps ) {
		
		delete next;
		next = cuda_jacobi_iteration( A, f, prev );
		norm = cuda_get_norm (next, prev);
		(*prev) = (*next);
		
		#ifdef CUDA_VERBOSE
		printf( "norm = %f\n", norm );
		#endif
	}
	
	return prev;
}

//==============================================================
//============ Single iteration ================================
vector * cuda_jacobi_iteration ( csr_matrix * A, vector * f, vector * prev ) {
#ifdef CUDA_VERBOSE
	printf("CUDA iteration started\n");
#endif

	int size = f->N;
	
	//cudaDeviceProp devProp;
	//cudaGetDeviceProperties ( &devProp, 0 );
	//int maxThreadsPerBlock = devProp.maxThreadsPerBlock;
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
	//printf("Pack memory\n");
	float * hostMatrixA = packMatrixA(A);
	float * hostVectorF = packVector(f);
	float * hostVectorPrev = packVector(prev);
	float * hostVectorNext = new float[size];
	//printf("End pack\n");
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
	
	vector * next = unpackVector ( hostVectorNext, size );
	
#ifdef CUDA_DEBUG
	for ( int i = 0; i < size; i++ )
		printf("next[%d]: %f\n", i, (*next)[i]);
#endif
	
	delete[] hostMatrixA;
	delete[] hostVectorF;
	delete[] hostVectorPrev;
	delete[] hostVectorNext;

	cudaFree(deviceMatrixA);
	cudaFree(deviceVectorF);
	cudaFree(deviceVectorPrev);
	cudaFree(deviceVectorNext);
	
#ifdef CUDA_VERBOSE
	printf("CUDA iteration ended\n");
#endif

	return next;
}

//==============================================================
//================ Common functions ============================
//==============================================================

//================ Norm of residual ============================
double cuda_get_norm ( vector * next, vector * prev ) {
	// —читаем норму, как максимальную поэлементную разность
	double result;
	result = fabs ( (*prev)[0] - (*next)[0] );
	
	for ( int i = 0; i < next->N; i++ ) {
		if ( fabs ( (*prev)[i] - (*next)[i] ) > result )
			result = fabs ( (*prev)[i] - (*next)[i] );
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

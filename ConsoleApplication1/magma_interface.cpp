#include "magma_interface.h"

vector * magma_jacobi_solve ( csr_matrix * A, vector * V ) {
	int size = V->N;
	int info;
	int * IPIV = new int[size];
	float * pA = packMatrixToBLAS(A);
	float *  pV = packVectorToBLAS(V);

	if ( magma_sgesv( size, 1, pA, size, IPIV, pV, size, &info ) != 0 )
		printf("Error in magma_sgesv\n");

	vector * R = unpackVectorFromBLAS(pV, size);

	deleteBlasMatrix(pA);
	deleteBlasVector(pV);
	delete[] IPIV;

	return R;
}

vector * magma_jacobi_solve ( float * M, float * V, int size ) {
	int * IPIV = new int[size];
	int info;
	
	printf("Calculating...\n");
	if ( magma_sgesv( size, 1, M, size, IPIV, V, size, &info ) != 0 )
		printf("Error in magma sgesv\n");
	
	printf("finished\nCopying to vector...\n");
	vector * R = unpackVectorFromBLAS( V, size );
	printf("done\n");
	delete[] IPIV;

	return R;
}

float * packMatrixToBLAS ( csr_matrix * input ) {
	int row = input->getRowCount();
	int col = input->getColCount();
	float * R = new float[row*col];
	
	for ( int i = 0; i < row; i++ )
		for ( int j = 0; j < col; j++ )
			R[i*col + j] = (float)input->getElement(i,j);
	
	return R;
}

float * packVectorToBLAS ( vector * input ) {
	float * R = new float[input->N];

	for ( int i = 0; i < input->N; i++ )
		R[i] = (float)((*input)[i]);

	return R;
}

vector * unpackVectorFromBLAS ( float * input, int size ) {
	vector * R = new vector(size);

	for ( int i = 0; i < size; i++ )
		(*R)[i] = (double)input[i];

	return R;
}

void deleteBlasVector ( float * input ) {
	delete[] input;
	input = NULL;
}

void deleteBlasMatrix ( float * R ) {
	delete[] R;
	R = NULL;
}
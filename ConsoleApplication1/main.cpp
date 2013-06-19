#include "main.h"


//=====================================================

int main ( int argc, char * argv[] ) {
//=====================================================

//================= Инициализация =====================
#ifdef MPI
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &PROC_NUMBER );
	MPI_Comm_rank( MPI_COMM_WORLD, &CURRENT_PROCESS );
#endif
#ifdef TEST
	#ifndef GENERATE_MATRIX
		#define GENERATE_MATRIX
	#endif
	#ifndef GENERATE_VECTOR
		#define GENERATE_VECTOR
	#endif
	#ifdef PRINT_RESULT
		#undef PRINT_RESULT
	#endif
	
	for ( DIMENSION = 2; DIMENSION < MAX_TEST_DIMENSION; DIMENSION *= DIMENSION_MULTIPLIER ) {
#elif DEFAULT_DIMENSION > 20
	#ifdef PRINT_RESULT
		#undef PRINT_RESULT
	#endif
#endif

#ifdef MPI
	if ( IS_ROOT_PROCESS ) {
#endif

		A = loadMatrixA();
		b = loadVectorB();
#ifdef MPI
	}

	sendMatrixToAllProcesses(&A);
	sendVectorToAllProcesses(&b);

#elif defined(MAGMA)
		printf("Packing matrix...\n");
		float * matrixForMagma = packMatrixToBLAS(A);
		printf("done\nPacking vector...\n");
		float * vectorForMagma = packVectorToBLAS(b);
		printf("done\n");
#endif
//=====================================================
//============== Решение системы ======================
//=====================================================

	timeStart = GetTickCount()/(double)1000;

#ifdef MPI
	if ( IS_SUCCESSIVE )
		result = jacobi_solve( A, b, JACOBI_EPSILON );
	else
		result = mpi_jacobi_solve( A, b, JACOBI_EPSILON, PROC_NUMBER, CURRENT_PROCESS );
#elif defined(CUDA)
		result = cuda_jacobi_solve( A, b, JACOBI_EPSILON );
#elif defined(MAGMA)
		result = magma_jacobi_solve( matrixForMagma, vectorForMagma, DIMENSION );
#else
		result = jacobi_solve( A, b, JACOBI_EPSILON );
#endif

	timeEnd = GetTickCount()/(double)1000;

#ifdef MPI
	if ( IS_ROOT_PROCESS )
#endif
		printTimeToFile("Time.txt");

#ifdef PRINT_RESULT
	check = (*A) * (*result);
#endif
	
//=====================================================
//================= Результаты ========================
//=====================================================
	
#ifdef MPI
	if ( IS_ROOT_PROCESS )
#endif
		printResult();

//=====================================================
//=============== Очищение памяти =====================
//=====================================================

#ifdef MAGMA
	deleteBlasMatrix(matrixForMagma );
	deleteBlasVector(vectorForMagma);
#endif
	
	memoryRelease();

#ifdef TEST
}
#endif

#ifdef MPI
	MPI_Finalize();
#endif

	return 0;
}

//==============================================================================================
//========== Вспомогательные функции ===========================================================
//==============================================================================================

void printResult ( void ) {

#ifdef PRINT_RESULT
	printf("===================\nInput Matrix\n");
	if (A) A->print();
	printf("===================\nInput F vector\n");
	if (b) b->print();
	printf("===================\nSolution\n");
	if (result) result->print();
	printf("===================\nMultiply\n");
	if (check) check->print();
#endif
	printf("===================\nTime\n");
	printf("time for dimension %d: %f\n", DIMENSION, timeEnd - timeStart);
	
	printf("===================\n");

}

//=====================================================

csr_matrix * loadMatrixA ( void ) {

#ifdef GENERATE_MATRIX
	printf("generating matrix...");
	sym_matrix * sym = new sym_matrix(DIMENSION);
	printf("random...");
	(*sym) = random_sym_matrix ( DIMENSION, MAX_MODULE_NUMBER, SINGULARITY_NUMBER);
	printf("copying to csr...");
	A = new csr_matrix(*sym);
	printf("record to file...\n");
	A->print("matrix.txt");
	delete sym;
	printf("matrix generated successfuly\n");
#else
	A = read_csr_matrix("matrix.txt");
#endif

	return A;
}

//=====================================================

vector * loadVectorB ( void ) {

#ifdef GENERATE_VECTOR	
	printf("generating vector...\n");
	b = random_vector(DIMENSION, MAX_MODULE_NUMBER);
	b->print("vector.txt");
	printf("vector generated successfuly\n");
#else
	b = read_vector("vector.txt");
#endif

	return b;
}

//=====================================================

void memoryRelease ( void ) {
	if (b) delete b;
	if (A) delete A;
	if (check) delete check;
	if (result) delete result;	
}

//=====================================================
//=========== Запаковка матрицы в массив ==============

double * packMatrix ( csr_matrix * input ) {
	int AAsize = input->getSizeOfAA();
	int JAsize = input->getSizeOfJA();
	int IAsize = (int)input->getSizeOfIA();
	int size = AAsize + JAsize + IAsize + 5;

	double * R = new double[size];
	int i;
	
	R[0] = (double)AAsize;
	R[1] = (double)JAsize;
	R[2] = (double)IAsize;
	R[3] = (double)input->getRowCount();
	R[4] = (double)input->getColCount();
	
	i = 5;
	for ( int k = 0; k < AAsize; k++, i++ )
		R[i] = input->getAAElement(k);
	
	for ( int k = 0; k < JAsize; k++, i++ )
		R[i] = (double)input->getJAElement(k);
	
	for ( int k = 0; k < IAsize; k++, i++ )
		R[i] = (double)input->getIAElement(k);

	return R;
}

//========== Распаковка матрицы из массива ==============

csr_matrix * unPackMatrix ( double * input ) {
	int AAsize =   (int)input[0];
	int JAsize =   (int)input[1];
	int IAsize =   (int)input[2];
	int rowCount = (int)input[3];
	int colCount = (int)input[4];

	csr_matrix * R = new csr_matrix( rowCount, colCount );
	
	int i = 5;
	for ( int k = 0; k < AAsize; k++, i++ )
		R->addAAElement( input[i] );
	
	for ( int k = 0; k < JAsize; k++, i++ )
		R->addJAElement( (int)input[i] );

	for ( int k = 0; k < IAsize; k++, i++ )
		R->addIAElement( (int)input[i] );
	
	return R;
}

//=====================================================
// Получение размера запакованного массива

int bufSize ( double * array, csr_matrix * matrix ) {
	return (matrix->getSizeOfAA() + matrix->getSizeOfJA() + matrix->getSizeOfIA() + 5);
}

//=====================================================
// Отсылка всем загруженной матрицы

void sendMatrixToAllProcesses( csr_matrix ** input ) {
	double * buf = NULL;
	int size;
	
	if ( IS_ROOT_PROCESS ) {
		buf = packMatrix(*input);
		size = bufSize(buf, *input);
	}
	
	MPI_Bcast(&size, 1, MPI_INT, ROOT_PROCESS, MPI_COMM_WORLD); // Рассылаем всем размер пересылаемого массива
	
	if ( IS_WORK_PROCESS )
		buf = new double[size];
	
	MPI_Bcast(buf,size, MPI_DOUBLE, ROOT_PROCESS, MPI_COMM_WORLD);

	if ( IS_WORK_PROCESS )
		*input = unPackMatrix(buf);

	//input->print();
}

//=====================================================
// Отсылка всем загруженного вектора
void sendVectorToAllProcesses( vector ** input ) {
	int size;

	if ( IS_ROOT_PROCESS )
		size = (*input)->N;

	MPI_Bcast(&size,1,MPI_INT, ROOT_PROCESS, MPI_COMM_WORLD); // Рассылка размера вектора всем процессам
	
	if ( IS_WORK_PROCESS )
		*input = new vector(size);

	MPI_Bcast((*input)->P,size,MPI_DOUBLE, ROOT_PROCESS, MPI_COMM_WORLD);
}

//=====================================================
//=====================================================
// Запись в файл времени выполнения
void printTimeToFile( char * filename ) {
	std::ofstream file;
	file.open("time.txt",std::ios::app);
	file << timeEnd - timeStart << std::endl;
	file.close();
}

//=====================================================
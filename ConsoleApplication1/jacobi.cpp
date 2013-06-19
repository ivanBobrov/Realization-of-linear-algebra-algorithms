
#include "jacobi.h"

//==============================================================
//=============== Якоби для формата симметричной матрицы =======
//==============================================================

vector * jacobi_solve ( sym_matrix * A, vector * f, double eps) {
	vector * prev = new vector(*f);
	vector next( (*f).N );
	double norm = eps;

	while ( norm >= eps ) {
		
		jacobi_iteration( A, f, prev, &next );
		norm = get_norm (&next, prev);
		(*prev) = next;
		
		#ifdef VERBOSE
		printf( "%f\n", norm );
		#endif
	}

	return prev;
}

vector * jacobi_iteration ( sym_matrix * A, vector * f, vector * prev, vector * next ) {
	// Просто вычисляем по формуле следующую итерацию
	for ( int i = 0; i < f->N; i++ ) {
		(*next)[i] = -(*f)[i];
		
		for ( int j = 0; j < f->N; j++ )
			if ( i != j )
				(*next)[i] += (*A).element(i,j) * (*prev)[j];
		
		(*next)[i] /= - (*A).element(i,i);
	}

	return next;
}

//==============================================================
//=============== Якоби для формата CSR матрицы ================
//==============================================================

vector * jacobi_solve ( csr_matrix * A, vector * f, double eps) {
	vector * prev = new vector(*f);
	vector next( (*f).N );
	double norm = eps;

	while ( norm >= eps ) {
		
		jacobi_iteration( A, f, prev, &next );
		norm = get_norm (&next, prev);
		(*prev) = next;
		
		#ifdef VERBOSE
		printf( "%f\n", norm );
		#endif
	}

	return prev;
}

vector * jacobi_iteration ( csr_matrix * A, vector * f, vector * prev, vector * next ) {
	// Просто вычисляем по формуле следующую итерацию
	for ( int i = 0; i < f->N; i++ ) {
		(*next)[i] = -(*f)[i];
		
		for ( int j = 0; j < f->N; j++ )
			if ( i != j )
				(*next)[i] += (*A).getElement(i,j) * (*prev)[j];
		
		(*next)[i] /= - (*A).getElement(i,i);
	}

	return next;
}

//==============================================================
//=============== MPI операции Якоби ===========================
//==============================================================

vector * mpi_jacobi_solve ( csr_matrix * A, vector * f, double eps, int PROC_NUMBER, int CURRENT_PROCESS ) {
	vector * prev = new vector(*f);
	vector * next = new vector( f->N );
	double norm = eps;

	while ( norm >= eps ) {
		
		mpi_jacobi_iteration( A, f, prev, next, PROC_NUMBER, CURRENT_PROCESS );
		norm = get_norm (next, prev);
		(*prev) = (*next);
		
		#ifdef VERBOSE
		printf( "%f\n", norm );
		#endif
	}

	return prev;
}

//==============================================================
// Итерация Якоби и распределение задач
vector * mpi_jacobi_iteration ( csr_matrix * A, vector * f, vector * prev, vector * next, int PROC_NUMBER, int CURRENT_PROCESS ) {

	double receivedData;
	int currentComponent = 0;
	int continueCalcTag = 1;
	MPI_Status status;

// ======= ROOT PROCESS ===========

	if ( IS_ROOT_PROCESS ) {
		
		while ( continueCalcTag ) {
			if ( MPI_SUCCESS != MPI_Recv(&receivedData, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status) )
				break;
		
			if ( status.MPI_TAG != next->N ) {
				(*next)[status.MPI_TAG] = receivedData;
				#ifdef MPI_VERBOSE
				printf("HOST received x[%d] = %f\n", status.MPI_TAG, receivedData);
				#endif
			}	
				
			if ( currentComponent < next->N ) {	
				MPI_Send( &currentComponent, 1, MPI_INT, status.MPI_SOURCE, 2, MPI_COMM_WORLD );
			}
			else if ( currentComponent >= next->N + PROC_NUMBER - 2 ) { // -1 от next->N: 0...next->N-1
																		// вторая -1 - число рабосих процессов меньше на управляющий
				continueCalcTag = 0;
			}
			currentComponent++;
		}
		
		for ( int i = 1; i < PROC_NUMBER; i++ )
			MPI_Send( &currentComponent, 1, MPI_INT, i, 0, MPI_COMM_WORLD );
		
	}

// ======= WORK PROCESS ===========

	if ( IS_WORK_PROCESS ) {
		
		receivedData = 0;
		MPI_Send ( &receivedData, 1, MPI_DOUBLE, ROOT_PROCESS, next->N, MPI_COMM_WORLD );
		
		while ( continueCalcTag ) {
			MPI_Recv ( &currentComponent, 1, MPI_INT, ROOT_PROCESS, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			if ( status.MPI_TAG == 0 ) {
				continueCalcTag = 0;
				break;
			}
			
			#ifdef MPI_VERBOSE
			printf("PROCESS %d received task %d\n", CURRENT_PROCESS, currentComponent );
			#endif
			// === непосредственно вычисления ===
			
			receivedData = calculateNextComponent( A, f, prev, currentComponent);
			
			// === ========================== ===
			MPI_Send ( &receivedData, 1, MPI_DOUBLE, ROOT_PROCESS, currentComponent, MPI_COMM_WORLD );
		}

	}

// ======= ALL PROCESS ===========	
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(next->P,next->N,MPI_DOUBLE,ROOT_PROCESS,MPI_COMM_WORLD);
	return next;
}

//==============================================================
// Вычисление одной компоненты итерации Якоби

double calculateNextComponent ( csr_matrix * A, vector * f, vector * prev, int component ) {
	double next;
	next = -(*f)[component];
		
	for ( int j = 0; j < f->N; j++ )
		if ( component != j )
			next += (*A).getElement(component,j) * (*prev)[j];
		
	next /= - (*A).getElement(component,component);
	
	return next;	
}

//==============================================================
//====== Общие функции последовательного Якоби =================
//==============================================================

double get_norm ( vector * next, vector * prev ) {
	// Считаем норму, как максимальную поэлементную разность
	double result;
	result = fabs ( (*prev)[0] - (*next)[0] );

	for ( int i = 0; i < next->N; i++ ) {
		if ( fabs ( (*prev)[i] - (*next)[i] ) > result )
			result = fabs ( (*prev)[i] - (*next)[i] );
	}

	return result;
}

//==============================================================

#ifdef GCC
// Если не определён, то определим модуль
double fabs ( double a ) {
	return (a>0) ? a : -a;
}
#endif

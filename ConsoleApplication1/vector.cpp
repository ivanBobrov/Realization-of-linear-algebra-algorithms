#include "vector.h"


//================================================================================================//
//=========================== Базовые функции класса вектор ======================================//
//================================================================================================//

vector::vector ( int num ) {
	if ( num < 1 ) {
		printf("warning: in vector constructor queried dimension is %d. Changed to 1\n",num);
		num = 1;		// Вектор должен быть положительной длины
	}
	
	N = num;			// Запись размера,
	P = new double[N];		// выделение памяти
	fill_zero();			// и заполнение нулями
	//printf("Create %x\n",this);
}

//==================================================================================

vector::vector ( const vector& init ) {
	N = init.N;
	P = new double[N];
	fill_zero();
	//printf("Create %x\n",this);
	for (int i = 0; i < N; i++ )
		P[i] = init.P[i];
}

//==================================================================================

vector::~vector ( void ) {
	//printf("Delete %x\n",this);
	delete P;		// Освобждение занятой памяти
}

//==================================================================================

void vector::fill_zero ( void ) {
	for (int i = 0; i < N; i++ )
		P[i] = 0;		// Заполнение нулями вектора
}

//================================================================================================//
//=========================== Перегрузка операторов класса вектор ================================//
//================================================================================================//

vector& vector::operator = ( vector const & right ) {
	if ( N != right.N ) {		// В случае различных размерностей
		delete P;		// удаление старой памяти и выделение новой
		N = right.N;		// с заполнением её нулями
		P = new double[N];	
		fill_zero();
	}

	for ( int i = 0; i < N; i++ )
		P[i] = right.P[i];	// Непосредственно присваивание
	return *this;
}

//==================================================================================

double& vector::operator[] ( int i ) {
	if ( i >= N ) {
		printf("Error: in vector operator[]. Index over of bounds. Queried %d while N = %d\n",i,N);
		return P[0];		// Проверка индексов
	} else
		return P[i];		// Возвращение ссылки на необходимый элемент
}

//==================================================================================

vector& vector::operator+= ( vector const & right ) {
	if ( N != right.N ) {
		printf("Error: in vector::operator+= : different dimensions\n");
		return *this;
	}
	
	for ( int i = 0; i < N; i++ )
		P[i] += right.P[i];
	
	return *this;
}

//==================================================================================

vector& vector::operator-= ( vector const & right ) {
	if ( N != right.N ) {
		printf("Error: in vector::operator-= : different dimensions\n");
		return *this;
	}
	
	for ( int i = 0; i < N; i++ )
		P[i] -= right.P[i];
	
	return *this;
}

//==================================================================================

vector operator + ( vector const & left, vector const & right ) {
	vector result(left.N);		// Вектор для записи результата
	
	if ( left.N != right.N ) {
		printf("Error: in vector::operator+ : different dimensions\n");
		return result; 		// Возврат нулевого вектора в случае ошибки в размерностях
	}
	
	for ( int i = 0; i < left.N; i++ )
		result.P[i] = left.P[i] + right.P[i];
					// Поэлементное сложение
	
	return result;
}

//==================================================================================

vector operator - ( vector const & left, vector const & right ) {
	vector result(left.N);		// Вектор для записи результата
	
	if ( left.N != right.N ) {
		printf("Error: in vector::operator- : different dimensions\n");
		return result; 		// Возврат нулевого вектора в случае ошибки в размерностях
	}
	
	for ( int i = 0; i < left.N; i++ )
		result.P[i] = left.P[i] - right.P[i];
					// Поэлементное сложение
	
	return result;
}

//==================================================================================

//================================================================================================//
//=========================== Общие функции класса вектор ================================//
//================================================================================================//

void vector::print ( const char * filename ) {
	if ( !strcmp(filename,"std") ) {
		printf("\n");
		for ( int i = 0; i< N; i++ )
			printf( " %7.3f ",P[i] );	
		printf("\n");
		
		return;
	}

	std::ofstream file(filename);
	file << N << std::endl;
	
	for (int i = 0; i < N; i++)
		file << P[i] << " ";
	
	file.close();

}

//==================================================================================

vector * read_vector ( const char * filename ) {
	vector * R;
	int dim;
	if ( !strcmp(filename, "std") ) {
		std::cout << "Vector input:\nInput length:" << std::endl;
		std::cin >> dim;
		std::cout << std::endl;
		R = new vector(dim);
		
		for ( int i = 0; i < dim; i++ ) {
			std::cout << "v[" << i << "]: ";
			std::cin >> R->P[i];
			std::cout << std::endl;
		}
		std::cout << std::endl;
	
	return R;
	}

	std::ifstream file(filename);
	
	if (!file.is_open()) {
		std::cout << "Error: in reading file (vector): no such file";
		return 0;
	}
	
	file >> dim;
	R = new vector(dim);
	
	for ( int i = 0; i < dim; i++ ) {
		file >> R->P[i];
	}

	file.close();

	return R;
}

//==================================================================================

vector * random_vector ( int dimension, int max ) {
	if ( dimension < 1 ) dimension = 1;
	vector * R = new vector(dimension);

	srand( time(NULL) );
	for ( int i = 0; i < dimension; i++ )
		R->P[i] = rand() % (2*max) - max;

	return R;
}
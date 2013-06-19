
#include "csr_matrix.h"


//=============== Конструкторы =====================
csr_matrix::csr_matrix ( int rows, int columns ) {
	AA.push_back(0);	// Добавляем первый нулевой элемент
	JA.push_back(0);	// для корректного отображения
	IA.push_back(0);
	
	row = rows;
	if ( columns == -1 )
		col = rows;
	else
		col = columns;
}

//====== Сжатие матрицы в формат CSR ======

csr_matrix::csr_matrix ( sym_matrix & A ) {
	row = A.dim;
	col = A.dim;
	int tmp = 0;		// Флаг. Если останется нулём, то исходная матрица нулевая.
	int IA_flag = -1;	// Флаг смены ряда. Для записи IA.
	
	for ( int i = 0; i < col; i++ ) {			// Бежим по всем
		for ( int j = 0; j < row; j++ ) {		// элементам
			
			if ( A.element(i,j) ) {							// Если элемент отличен от нуля
				AA.push_back( A.element(i,j) );				// добавляем его значение в AA
				JA.push_back(j);							// добавляем номер строки в JA
				if ( IA_flag < i ) { 						// Если сменился ряд, то
					IA_flag = i;							// пишем индекс элемента
					IA.push_back(AA.size());				// в IA.
				}
				
				if (!tmp) tmp = 1;						// Хотя бы один элемент ненулевой.
			}
			
		}
	}
	
	if (tmp) IA.push_back( AA.size() + 1 );				// Если исходная матрица ненулевая,
		else {											// добавляем в IA завершающее
			AA.push_back(0);							// значение AA
			JA.push_back(0);							// (т.к. мы сделали clear для всех
			IA.push_back(0);							// этих массивов)
		}	
}

//============== Private методы =====================

void csr_matrix::clearArrays ( void ) {
	AA.clear();
	IA.clear();
	JA.clear();
}

//============== Public методы =====================

//======== Взятие элемента ========
double csr_matrix::getElement ( int i, int j ) {
	if ( (i > row) || (i < 0) || (j > col) || (j < 0) ) {
		printf("Error: Trying to get element outside matrix's bounds (CSR)\n");
		return 0;
							// Если запрошенный элемент превышает размерность
	}						// то возвращаем ноль.
	
	for ( int k = (IA[i] - 1); k < (IA[i+1] - 1); k++ ) 
		if ( JA[ k ] == j )				// Проходим по всей i-ой строке,
			return AA[ k ];				// сравнивая совпадает ли у элемента
										// номер столбца (j) и JA.
	return 0;					
	
	printf("Error: nothing to return in csr_matrix element");
}

//=================================

//=========== Вывод матрицы =======
void csr_matrix::print ( const char * filename ) {
	if ( !strcmp (filename,"std") ) {
		printf("\n");
		for ( int i = 0; i < row; i++ ) {
			for ( int j = 0; j < col; j++ )
				printf(" %7.3f ",this->getElement(i,j));
			printf("\n");
		}
		printf("\n");
		return;
	}
	

	std::ofstream file(filename);
	file << row << " " << col << " " << AA.size() << " " << JA.size() << " " << IA.size() << std::endl;
	
	for (int i = 0; i < (int)AA.size(); i++)
		file << AA[i] << " ";
	file << std::endl;
	for (int i = 0; i < (int)JA.size(); i++)
		file << JA[i] << " ";
	file << std::endl;
	for (int i = 0; i < (int)IA.size(); i++)
		file << IA[i] << " ";
	
	file.close();
}

//===================================

//=========== Ввод матрицы =========

csr_matrix * read_csr_matrix ( const char * filename ) {
	int col, row;
	csr_matrix * R;

	if ( !strcmp (filename,"std") ) {
		
		std::cout << "Col number: "   << std::endl;
		std::cin  >> col;
		std::cout << "Row number: "   << std::endl;
		std::cin  >> row;
		std::cout << "Input matrix: " << std::endl;
	
		R = new csr_matrix(row, col);
		R->clearArrays();

		double temp_input;
		bool new_line = true;
		for ( int i = 0; i < row; i++ ) {
			
			for ( int j = 0; j < col; j++ ) {
				std::cin >> temp_input;
				
				if (temp_input) {
					R->AA.push_back(temp_input);
					R->JA.push_back(j);
				}
				if (new_line) {
					new_line = false;
					R->IA.push_back(R->AA.size());
				}
			}		
			
			new_line = true;
		}
		R->IA.push_back(R->AA.size() + 1); // Завершающее значение. Чтобы циклы в print обрывался.

	return R;
	}

	std::ifstream file(filename);
	
	if (!file.is_open()) {
		std::cout << "Error: in reading file: no such file";
		return 0;
	}
	
	int AA_size, JA_size, IA_size;
	file >> row >> col >> AA_size >> JA_size >> IA_size;
	R = new csr_matrix(row, col);
	R->clearArrays();

	double temp_value;
	int temp_int_value;
	for ( int i = 0; i < AA_size; i++ ) {
		file >> temp_value;
		R->AA.push_back(temp_value);
	}
	for ( int i = 0; i < JA_size; i++ ) {
		file >> temp_int_value;
		R->JA.push_back(temp_int_value);
	}
	for ( int i = 0; i < IA_size+1; i++ ) {
		file >> temp_int_value;
		R->IA.push_back(temp_int_value);
	}

	file.close();
	return R;
}

//=========== Общие функции =========

vector* operator * ( csr_matrix & matrix, vector const & vect ) {
	int row = matrix.getRowCount();
	int col = matrix.getColCount();

	vector * R = new vector(row);
		if ( col != vect.N ) {
		printf("Error: in csr_matrix::operator*(to vector) : different dimensions [%d,%d]*[%d]\n", row, col, vect.N);
		return R;
	}
	
	for ( int i = 0; i < row; i++ )
		for ( int j = 0; j < col; j++ )
			R->P[i] += matrix.getElement(i,j) * vect.P[j];

	return R;
}

//===================================
int csr_matrix::getRowCount( void ) { return row; }
int csr_matrix::getColCount( void ) { return col; }
int csr_matrix::getSizeOfAA( void ) { return AA.size(); }
int csr_matrix::getSizeOfJA( void ) { return JA.size(); }
int csr_matrix::getSizeOfIA( void ) { return IA.size(); }
int csr_matrix::getJAElement( int i ) { return JA[i]; }
int csr_matrix::getIAElement( int i ) { return IA[i]; }
double csr_matrix::getAAElement( int i ) { return AA[i]; }
void csr_matrix::setAAElement( int place, double value ) { AA[place] = value; }
void csr_matrix::setJAElement( int place, int value ) { JA[place] = value; }
void csr_matrix::setIAElement( int place, int value ) { IA[place] = value; }

void csr_matrix::addAAElement( double value ) {
	if( (AA.size() == 1) && (AA[0] == 0) )
		clearArrays();	
	AA.push_back(value);
}

void csr_matrix::addJAElement( int value ) { JA.push_back(value); }
void csr_matrix::addIAElement( int value ) { IA.push_back(value); }

//===================================

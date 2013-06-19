#include "sym_matrix.h"

//================================================================================================//
//=========================== Базовые функции класса симметричных матриц =========================//
//================================================================================================//

sym_matrix::sym_matrix ( int N = 1 ) {
	if (N < 1) N = 1;
	dim = N;
	len = N*(N+1)/2;
	A = new double[len];
	for (int i = 0; i < len; i++ )
		A[i] = 0;
}

//==================================================================================

sym_matrix::sym_matrix ( sym_matrix const & init ) {
	dim = init.dim;
	len = init.len;
	A = new double[len];
	for (int i = 0; i < len; i++ )
		A[i] = init.A[i];
}

//==================================================================================

sym_matrix::~sym_matrix ( void ) {
	delete A;
}

//================================================================================================//
//=========================== Перегрузка операторов класса симметричных матриц ===================//
//================================================================================================//

//============================//
//======= Члены класса =======//=======================================
//============================//

sym_matrix& sym_matrix::operator = ( sym_matrix const & right ) {
	if ( dim != right.dim ) {
		dim = right.dim;
		len = right.len;
		delete A;
		A = new double[len];
	}
	
	for ( int i = 0; i < len; i++ )
		A[i] = right.A[i];
	
	return *this;
}

//==================================================================================

double& sym_matrix::element ( int row, int col ) {
	if ( col >= row )
		return A[ ((col+1)*col/2) + row ];
	else
		return A[ ((row+1)*row/2) + col ];
}

//==================================================================================

sym_matrix& sym_matrix::operator += ( sym_matrix const & right ) {
	if ( dim != right.dim ) {
		printf("Error: in sym_matrix::operator+= : different dimensions\n");
		return *this;
	}
	
	for ( int i = 0; i < len; i++ )
		A[i] +=right.A[i];
	
	return *this;
}

//==================================================================================

sym_matrix& sym_matrix::operator -= ( sym_matrix const & right ) { return *this += -right; }

//==================================================================================

sym_matrix sym_matrix::operator- ( void ) const {
	sym_matrix result(dim);
	
	for ( int i = 0; i < len; i++ )
		result.A[i] = -A[i];
	
	return result;
}

//==================================================================================

sym_matrix& sym_matrix::operator *= ( sym_matrix & right ) {
	sym_matrix result(dim);
	
	if ( dim != right.dim ) {
		printf("Error: in sym_matrix::operator*= : different dimensions\n");
		return *this;
	}
	
	for ( int i = 0; i < dim; i++ )
		for ( int j = i; j < dim; j++ )
			for ( int k = 0; k < dim; k++ )
				result.element(i,j) += this->element(i,k)*right.element(k,j);
			
	*this = result;
	
	return *this;
}

//============================//
//==== Функции вне класса ====//=======================================
//============================//

sym_matrix operator + ( sym_matrix const & left, sym_matrix const & right ) {
	if ( left.dim != right.dim ) {
		printf("Error: in sym_matrix::operator+ : different dimensions [%d]+[%d]\n",left.dim,right.dim);
		return 0;
	}
	
	sym_matrix result(left.dim);
	
	for ( int i = 0; i < left.len; i++ )
		result.A[i] = left.A[i] + right.A[i];
	
	return result;
}

//==================================================================================

sym_matrix operator + ( sym_matrix const & left, double number ) {
	sym_matrix result(left.dim);
	
	for (int i = 0; i < left.len; i++)
		result.A[i] = left.A[i] + number;
	
	return result;
}

//==================================================================================

sym_matrix operator + ( double number, sym_matrix const & right ) { return right + number; }

//==================================================================================
//==================================================================================

sym_matrix operator - ( sym_matrix const & left, sym_matrix const & right ) { return left + (-right); }

//==================================================================================

sym_matrix operator - ( sym_matrix const & left, double number ) { return left + (-number); }

//==================================================================================

sym_matrix operator - ( double number, sym_matrix const & right ) { return (-right) + number; }

//==================================================================================
//==================================================================================

sym_matrix operator * ( sym_matrix & left, sym_matrix & right ) {
	if ( left.dim != right.dim ) {
		printf("Error: in sym_matrix::operator* : different dimensions\n");
		return 0;
	}
	
	sym_matrix result(left.dim);
	
	for ( int i = 0; i < left.dim; i++ )
		for ( int j = i; j < left.dim; j++ )
			for ( int k = 0; k < left.dim; k++ )
				result.element(i,j) += left.element(i,k)*right.element(k,j);
	
	return result;
}

//==================================================================================

sym_matrix operator * ( sym_matrix const & left, double number ) {
	sym_matrix result(left.dim);
	
	for (int i = 0; i < left.len; i++)
		result.A[i] = left.A[i] * number;
	
	return result;
}

//==================================================================================

sym_matrix operator * ( double number, sym_matrix const & right ) { return right*number; }

//==================================================================================

vector operator * ( sym_matrix & left, vector const & vect ) {
	vector result(left.dim);
	if ( left.dim != vect.N ) {
		printf("Error: in sym_matrix::operator*(to vector) : different dimensions [%d]*[%d]\n",left.dim,vect.N);
		return result;
	}
	
	for ( int i = 0; i < left.dim; i++ )
		for ( int j = 0; j < left.dim; j++ )
			result.P[i] += left.element(i,j)*vect.P[j];

	return result;

}

//==================================================================================

std::vector<double> operator * ( sym_matrix & left, std::vector<double> vect ) {
	std::vector<double> result(left.dim,0);
	if ( left.dim != vect.size() ) {
		printf("Error: in sym_matrix::operator*(to vector) : different dimensions [%d]*[%d]\n",left.dim,vect.size() );
		return result;
	}
	
	for ( int i = 0; i < left.dim; i++ )
		for ( int j = 0; j < left.dim; j++ )
			result[i] += left.element(i,j)*vect[j];
		
	return result;
}

//==================================================================================
//==================================================================================

sym_matrix operator / ( sym_matrix const & left, double number ) {
	sym_matrix result(left.dim);
	
	for (int i = 0; i < left.len; i++)
		result.A[i] = left.A[i] / number;
	
	return result;
}

//================================================================================================//
//=========================== Общие функции класса симметричных матриц ===========================//
//================================================================================================//

int sym_matrix::size( void ) { return dim; }

//==================================================================================

void sym_matrix::print ( const char* filename ) {
	if ( !strcmp (filename,"stdo") ) {
		printf("\n");
		for ( int i = 0; i < dim; i++ ) {
			for ( int j = 0; j < dim; j++ )
				printf(" %7.3f ",this->element(i,j));
			printf("\n");
		}
		printf("\n");
		return;
	}
	
	std::ofstream file(filename);
	file << dim << std::endl;
	/*
	for ( int i = 0; i < dim; i++ ) {
		for ( int j = 0; j < dim; j++ )
			file << this->element(i,j) << " ";
		file << std::endl;
	}*/
	
	for ( int i = 0; i < len; i++ )
			file << A[i] << " ";
	
	file.close();
}

//================================================================================================//
//=========================== Общие функции ======================================================//
//================================================================================================//

sym_matrix read_sym_matrix ( const char* filename ) {
	if ( !strcmp (filename,"stdo") ) {
		int N;
		std::cout << "Размерность: ";
		std::cin >> N;
		sym_matrix result(N);
		
		for (int i = 0; i < result.len; i++ )
			std::cin >> result.A[i];
		
		return result;
	}
	
	int N;
	std::ifstream file(filename);
	
	if (!file.is_open()) {
		std::cout << "Error: in reading file: no such file";
		return 0;
	}
	
	file >> N;
	sym_matrix result(N);
	for (int i = 0; i < result.len; i++ )
		file >> result.A[i];
	
	file.close();
	return result;
}

//==================================================================================

sym_matrix random_sym_matrix ( int dim, int max_num, double singul ) {
	if ( dim < 1 ) dim = 1;
	sym_matrix result(dim);

	srand( time(NULL) );
	for ( int i = 0; i < result.dim; i++)
		for ( int j = i; j < result.dim; j++ ) {
			if ( !(rand() % 50) )
				result.element(i,j) = rand() % (2*max_num) - max_num;
			else
				result.element(i,j) = 0;
			//result.element(i,j)*= normal_func(j-i,dim-1,(double)dim*1/singul);
		}
	
	for ( int i = 0; i < result.dim; i++ ) {
		result.element(i,i) = module (result.element(i,i) );
		for ( int j = 0; j < result.dim; j++ )
			result.element(i,i) += module( result.element(i,j) );
	if (!result.element(i,i))
		result.element(i,i) = rand() % (2*max_num) - max_num;
	}

	return result;
}

//==================================================================================

double normal_func ( double diff, double max_diff, double point ) {
	if ( diff > point )
		return (diff-max_diff)*point/(max_diff*(point-max_diff));
	else
		return ((point-max_diff)/(max_diff*point))*diff + 1;
}

//==================================================================================

double module ( double N ) { return N>0 ? N : -N;}
#ifndef SYM_MATRIX_H
#define SYM_MATRIX_H


#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <ctime>

#include "vector.h"

class sym_matrix {
	public:
	//private:
		double *A;		// Массив значений
		int dim,len;		// dim - размерность матрицы; len - длина массива данных;
	//public:
		sym_matrix ( int );
		~sym_matrix ( void );
		sym_matrix ( sym_matrix const & );
		
		double& element ( int row, int col );
		sym_matrix& operator  = ( sym_matrix const & );
		sym_matrix& operator += ( sym_matrix const & );
		sym_matrix& operator -= ( sym_matrix const & );
		sym_matrix& operator *= ( sym_matrix       & );
		sym_matrix  operator  - (     void     ) const;
		
		int size( void );
		void print ( const char* filename = "stdo" );

};

sym_matrix operator + ( sym_matrix const &, sym_matrix const & );
sym_matrix operator + ( sym_matrix const &, double );
sym_matrix operator + ( double, sym_matrix const & );

sym_matrix operator - ( sym_matrix const &, sym_matrix const & );
sym_matrix operator - ( sym_matrix const &, double );
sym_matrix operator - ( double, sym_matrix const & );

sym_matrix operator * ( sym_matrix       &, sym_matrix       & );
sym_matrix operator * ( sym_matrix const &, double );
sym_matrix operator * ( double, sym_matrix const & );

vector     operator * ( sym_matrix       &, vector const & );
std::vector<double> operator * ( sym_matrix &, std::vector<double> );

sym_matrix operator / ( sym_matrix const &, double );

sym_matrix read_sym_matrix ( const char* filename = "stdo" );
sym_matrix random_sym_matrix ( int dimension, int max_num, double singul = 1.0 );
double normal_func ( double diff, double max_diff, double point );

double module (double);
#endif

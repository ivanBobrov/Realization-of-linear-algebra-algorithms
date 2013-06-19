#ifndef CSR_MATRIX_H
#define CSR_MATRIX_H

#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "sym_matrix.h"

class csr_matrix {
	private:
		int col;
		int row;
		
		std::vector<double> AA;
		std::vector<int> JA;
		std::vector<int> IA;

		void clearArrays( void );
	
	public:	
		csr_matrix ( int rows, int columns = -1 );
		csr_matrix ( sym_matrix & A );
		
		double getElement ( int i, int j );
		void print ( const char * filename = "std" );

		int getRowCount( void );
		int getColCount( void );
		int getSizeOfAA( void );
		int getSizeOfJA( void );
		int getSizeOfIA( void );

		double getAAElement( int );
		int getJAElement( int );
		int getIAElement( int );
		void setAAElement( int place, double value );
		void setJAElement( int place, int value );
		void setIAElement( int place, int value );
		void addAAElement( double value );
		void addJAElement( int value );
		void addIAElement( int value ); 

		friend csr_matrix * read_csr_matrix ( const char * filename = "std" );

/*
 ласс будет некорректно работать в случае, если в матрице содержитс€ нулева€ строка.
–абота с такими матрицами не осуществл€етс€.
*/

};

vector * operator * ( csr_matrix &, vector const & );
#endif
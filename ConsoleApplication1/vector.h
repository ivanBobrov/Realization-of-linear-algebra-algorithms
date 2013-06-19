#ifndef VECTOR_H
#define VECTOR_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <ctime>

class vector {
	public:
	int N;
	double *P;
	
	vector ( int num );
	vector ( const vector& );
	~vector ( void );
	
	vector& operator = ( vector const & right );
	vector& operator+= ( vector const & right );
	vector& operator-= ( vector const & right );
	double& operator[] ( int );
	
	void fill_zero ( void ); // Заполнение вектора нулями
	void print ( const char * filename = "std" );
};

vector operator + ( vector const & left, vector const & right );
vector operator - ( vector const & left, vector const & right );
vector * read_vector ( const char * filename = "std" );
vector * random_vector ( int dimension, int max_number );
#endif
#include "headers.hpp"

//niveau blas1 omp_parallele
double BLAS1(Vector a, Vector b)
{
	double s = 0;
	int i;

	#pragma omp parallel
	{
		#pragma omp for reduction(+:s)
		for( i = 0; i < a.rows(); i++)
	  		s += a(i)*b(i);
	}
  return s;
}

Vector BLAS2(Matrix A, Vector x)
{
	Vector s = Vector::Zero(A.rows());
  	Vector a = Vector::Zero(x.rows());
		int i,j;

	  	for ( i = 0; i < A.rows(); ++i)
	  	{
	  		for (j = 0; j < x.rows(); ++j)
	  			a(j) = A(i,j);
	  			s(i) = BLAS1(a, x);
	  	}

  return s;
}

Vector BLAS3(Matrix A, Vector x)
{

	Vector s = BLAS2(A, x);

  s = BLAS2(A, s) + s + x;

  // s.conservativeResize(x.rows());

  return s;
}

#include "headers.hpp"


void PuissanceIt(Problem& p, double tol)
{
	p.v = p.v/sqrt(BLAS1(p.v,p.v));
	Vector y = BLAS2(p.A,p.v);
	Vector x = y/sqrt(BLAS1(y,y));
	double err = abs(abs(BLAS1(p.v,x)) - 1);
	while (err > tol)
	{
		p.v = x;
		y = BLAS2(p.A,p.v);
		x = y/sqrt(BLAS1(y,y));
		err = abs(abs(BLAS1(p.v,x)) - 1);
	}
	p.lambda = y(0)/p.v(0);
}

void Deflation(Problem& p, Matrix& B, int m, double tol)
{
	Vector v = p.v;
	B = Matrix::Zero(p.v.rows(), m);

	p.Lmbd = Vector::Zero(m);

	PuissanceIt(p, tol);
	p.Lmbd(0) = p.lambda;
	for (int i = 0; i < p.v.rows(); ++i)
		B(i,0) = p.v(i);

	for (int k = 1; k < m; ++k)
	{
		Matrix M = Prodvc(p.v, p.v);
		p.A = p.A - p.lambda*M;
		p.v = v;
		PuissanceIt(p, tol);
		p.Lmbd(k) = p.lambda;
		for (int i = 0; i < p.v.rows(); ++i)
			B(i,k) = p.v(i);
	}
}

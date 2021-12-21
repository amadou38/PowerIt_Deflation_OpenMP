#include "headers.hpp"



void Tests(Problem& p, int cas, int dim)
{
	if (cas == 1)
	{
		double a = 11111111, b = 90909090, c = 10891089, d = 8910891, e = 11108889, f = 9089091, g = 10888911, h = 8909109;

		p.A.resize(8,8);
		p.Lmbd.resize(8);
		p.v = Vector::Random(8);
		Matrix A = Matrix::Zero(2,2), B = A, C = A, D = A;

		A(0,0) = a; A(0,1) = -b; A(1,0) = -b; A(1,1) = a;
		B(0,0) = -c; B(0,1) = d; B(1,0) = d; B(1,1) = -c;
		C(0,0) = -e; C(0,1) = f; C(1,0) = f; C(1,1) = -c;
		D(0,0) = g; D(0,1) = -h; D(1,0) = -h; D(1,1) = g;
		Block(p.A, A, 0,2, 0,2); Block(p.A, B, 0,2, 2,4); Block(p.A, C, 0,2, 4,6); Block(p.A, D, 0,2, 6,8);
		Block(p.A, B, 2,4, 0,2); Block(p.A, A, 2,4, 2,4); Block(p.A, D, 2,4, 4,6); Block(p.A, C, 2,4, 6,8);
		Block(p.A, C, 4,6, 0,2); Block(p.A, D, 4,6, 2,4); Block(p.A, A, 4,6, 4,6); Block(p.A, B, 4,6, 6,8);
		Block(p.A, D, 6,8, 0,2); Block(p.A, C, 6,8, 2,4); Block(p.A, B, 6,8, 4,6); Block(p.A, A, 6,8, 6,8);

		for (int i = 0; i < p.A.rows(); ++i)
		    p.Lmbd(i) = 8*pow(10,i);
	}
	if (cas == 2)
	{
		p.A = Matrix::Zero(dim,dim);
		p.Lmbd = Vector::Zero(dim);
		p.v = Vector::Random(dim);

		for (int j = 0; j < dim; ++j)
			for (int i = 0; i <= j; ++i)
				p.A(i,j) = dim + 1 - j;
		for (int j = 0; j < dim-1; ++j)
			for (int i = 0; i <= j+1; ++i)
				p.A(i,j) = dim + 1 - i;

		for (int i = 0; i < p.A.rows(); ++i)
		    p.Lmbd(i) = 1/(2 - 2*cos((2*i-1)*M_PI/(2*dim+1)));
	}
	if (cas == 3)
	{
		p.A = Matrix::Zero(dim,dim);
		p.Lmbd = Vector::Zero(dim);
		p.v = Vector::Random(dim);

		double a = 10, b = 4;

		for (int i = 0; i < dim; ++i)
		{
			p.A(i,i) = a;
			if (i < dim-1)
				p.A(i,i+1) = b;
			if (i > 0)
				p.A(i,i-1) = b;
		}

		for (int i = 0; i < dim; ++i){
		    p.Lmbd(i) = a - 2*b*cos((i+1)*M_PI/(dim+1));
		    // p.v(i) = i+1;
		}
	}
	if (cas == 4)
	{
		p.A.resize(dim,dim);
		p.Lmbd.resize(dim);
		p.v = Vector::Random(dim);

		double a = 11111111, b = 90909090;

		for (int i = 0; i < dim; ++i)
		{
			p.A(i,i) = a;
			if (i < dim-1)
				p.A(i,i+1) = b;
			if (i > 0)
				p.A(i,i-1) = b;
		}

		for (int i = 0; i < dim; ++i)
		    p.Lmbd(i) = a - 2*b*cos((i+1)*M_PI/(dim+1));
	}

}

void sortVect(Vector& v)
{
	double max;

	for (int i = 0; i < v.rows()-1; ++i)
		for (int j = 0; j < v.rows()-1-i; ++j)
			if (v(j) < v(j+1))
			{
				max = v(j+1);
				v(j+1) = v(j);
				v(j) = max;
			}
}

void sort(Eigen::MatrixXd& Eigval, Eigen::MatrixXd& Eigvect)
{
	for (int i = 0; i < Eigval.cols()-1; ++i)
		for (int j = 0; j < Eigval.cols()-1-i; ++j)
			if (Eigval(j,j) < Eigval(j+1,j+1))
			{
				double max = Eigval(j+1,j+1);
				Eigval(j+1,j+1) = Eigval(j,j);
				Eigval(j,j) = max;
				for (int k = 0; k < Eigvect.rows(); ++k)
				{
					max = Eigvect(k,j+1);
					Eigvect(k,j+1) = Eigvect(k,j);
					Eigvect(k,j) = max;
				}
			}
}

void MatrixCols(Vector& v, Matrix& A, int j)
{
	for (int i = 0; i < A.cols(); ++i)
		v(i) = A(i,j);
}

void MatrixRows(Vector& v, Matrix& A, int i)
{
	for (int j = 0; j < A.rows(); ++j)
		v(j) = A(i,j);
}


Matrix Prodvc(Vector a, Vector b)
{
	Matrix S = Matrix::Zero(a.rows(), b.rows());

#pragma omp parallel
{
	#pragma for schedule(runtime)
	for(int j = 0; j < b.rows(); j++)
		for(int i = 0; i < a.rows(); ++i)
			S(i,j) = a(i)*b(j);
}

  	return S;
}

void Block(Matrix& AA, Matrix& B, int i1, int i2, int j1, int j2)
{
	int k = 0, l = 0;
	for (int i = i1; i < i2; ++i)
	{
		l = 0;
		for (int j = j1; j < j2; ++j)
		{
			AA(i,j) = B(k,l);
			l++;
		}
		k++;
	}
}

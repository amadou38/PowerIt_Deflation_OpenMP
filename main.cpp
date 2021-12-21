#include "headers.hpp"
#include <omp.h>

int main(int argc, char *argv[])
{
	cout.precision(2);

	Problem p;

    int n = 24;   // taille de la matrice

    Tests(p, 3, n);
    Matrix A = p.A;
    Vector v = p.v;

    int m = 6;
    double tol = 1e-15;
    Matrix B;
    Vector Lambda;

    sortVect(p.Lmbd);

    cout << "Data present in the file: " << endl << endl;
    cout << "Matrix p.A (global): \n" << A << endl << endl;
    cout << "Initial vector p.v (global): \n" << v << endl << endl;
    cout << "\nExact Eigenvalues : " << endl;
    for (int i = 0; i < p.Lmbd.rows(); i++)
        cout << p.Lmbd(i) << "  ";
    cout << endl;

    auto t1 = omp_get_wtime();
    // Deflation method (with Puissance It)
    Deflation(p, B, m, tol);
    auto t2 = omp_get_wtime();

    auto T = t2-t1;

    cout << "Approximate Eigenvalues (for Deflation): " << endl;
    for (int i = 0; i < m; i++)
        cout << p.Lmbd(i) << "  ";
    cout << "\nAnd associated eigenvectors (global):\n" << B << endl << endl;

    cout << "\n\nPuissance It & Deflation runtime: " << T << " sec\n\n";

	return 0;
}

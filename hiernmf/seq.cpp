//
// @author John Farrell
// @date 2/2/20
//

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;


// returns the maximum eigenvalue of A-transpose * A
vec powIter (mat A, vec q);


/*
 * This code runs the power iteration using the armadillo library in C++
 * The power iteration calulates the maximum eigenvalue in some matrix A
 * This code in particular will calculate the max(eig(A'*A))
 *
 * Later we will convert this into MPI-based parallelism
 */
int main () {

	// UNIT TESTS
	mat A = randu<mat>(5,5);
	vec q = randu<vec>(5);
	cout << A.n_rows << "x" << A.n_cols << " matrix A'A: " << endl;
	cout << A.t()*A << endl;
	cout << "Max singular value: \n" << max(svd(A.t()*A)) << endl;
	cout << (powIter(A,q) - max(svd(A.t()*A))) << endl;	
	return 0;
}


vec powIter (mat A, vec q) {
	
	// normalize q beforehand
	q = q / norm(q, 2);

	vec lambda, l2;
	vec z;
	int i = 0;
	double epsilon = 1.0;

	cout << "\nConvergence: " << endl;

	// powerIteration
	while (i < 10 && epsilon > 0.000001) {
		z = A * q;
		q = A.t() * z;
		q = q / norm(q, 2);

		// store previous lambda for convergence measuring
		if (i != 0) {
			l2 = lambda;
			lambda = z.t() * z;
			epsilon = as_scalar(lambda - l2);
		}
		else {
			lambda = z.t() * z;
		}
		
		// track convergence
		cout << "lambda:  " << lambda << endl;
		i++;
	}

	return lambda;
}


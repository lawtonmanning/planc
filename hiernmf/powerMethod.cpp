//
// @author John Farrell
// @date 2/2/20
//
// 2/7/20 - begin to parallelize the method
//

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;


// returns the maximum eigenvalue of A-transpose * A
double powIter (mat A);


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
	cout << A.n_rows << "x" << A.n_cols << " matrix A: " << endl;
	cout << A << endl;
	cout << "Max singular value: \n" << max(svd(A)) << endl;
	cout << (powIter(A) - max(svd(A))) << endl;	
	return 0;
}


/* NOTES:
 * 
 * 	Is there a way to change lambda to auto so that it does not need to convert from vector? Hard code as double?
 *	 
 */
double powIter (mat A) {
	
	// generate some random q vector
	vec q = randu<vec>(size(A,1));
	
	// normalize q beforehand
	q = q / norm(q, 2);

	vec z;
	auto lambda = norm(q,2);
	auto l2 = lambda;
	int i = 0;
	double epsilon = 1.0;

	cout << "\nConvergence: " << endl;

	// powerIteration
	while (i < 10 && epsilon > 0.00001) {
		z = A * q;
		q = A.t() * z;
		q = q / norm(q, 2);

		// store previous lambda for convergence measuring
		if (i == 0)
			lambda = norm(z,2);
		else {
			l2 = lambda;
			lambda = norm(z,2);
			epsilon = abs(lambda - l2)/(lambda);
		}
		
		// track convergence
		cout << "lambda:  " << lambda << endl;
		i++;
	}

	return lambda;
}


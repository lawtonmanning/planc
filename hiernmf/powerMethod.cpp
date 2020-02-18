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


/*
 *	Power Method : Parallel Implementation
 *		- This code receives some matrix A, and uses the Power Method to converge to the largest singular value
 *		- In this case we use A' * A to find the principal sigma value, squared
 *	
 */
double powIter (mat A) {
	
	// generate some random q vector and normalize
	vec q = randu<vec>(size(A,1));
	q = q / norm(q, 2);

	// prepare variables
	vec z;
	auto sigma = norm(q,2);
	auto s2 = sigma;
	double epsilon = 1.0;

	cout << "\nConvergence: " << endl;

	// converge to sigma1 using the power method
	int iter = 0;
	while (iter < 10 && epsilon > 0.00001) {
		z = A * q;	
		z = z / norm(z,2);	// normalize z to manage large numbers
		q = A.t() * z;		
		sigma = norm(q,2);	// set sigma here to prevent repetitive calculations
		q = q / sigma;
			
		epsilon = abs(sigma - s2)/(sigma);
		s2 = sigma;
		iter++;

		// debugging output
		cout << "Sigma:  " << sigma << endl;
	}

	return sigma;
}


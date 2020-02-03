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

	// Notes:
	// eig_gen(X) is the general eigenvalue decomp for dense matrices 
	// 	(non-sym | non-hermitian)

	// UNIT TESTS
	mat A = randu<mat>(5,5);
	vec q = randu<vec>(5);
	cout << A.t()*A << "\n" << svd(A.t()*A) << endl;
	cout << max(svd(A.t()*A)) << endl;
	cout << "Accuracy A = 5x5:\n" << (powIter(A, q) - max(svd(A.t()*A))) << endl;
	
	return 0;
}


// Notes:
// 	add verbose to track sequence?
vec powIter (mat A, vec q) {
	
	// normalize q
	q = q / norm(q, 2);

	vec lambda;
	int i = 0;

	// temporary -- change into discrete steps without explicitly multiplying A'A
	while (i < 10) {
		q = A.t() * A * q;
		q = q / norm(q, 2);
		lambda = q.t() * A.t() * A * q;
		
		// track convergence
		cout << "lambda:  " << lambda << endl;
		i++;
	}

	return lambda;
}


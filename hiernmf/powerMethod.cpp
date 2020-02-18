//
// @author John Farrell
// @date created 2/2/20
//
// 2/17/20
// Minor improvements to power method
// 	- Changed lambda to sigma
// 	- Removed unnecessary if/else statement
// 	- normalized vector z between calculations
// MPI additions
// 	- MPI library included
// 	- AllREDUCE placed after q is calculated the first time
//
// 2/../..
// SUMMARY
// 	- POINTS
//

#include <iostream>
#include <armadillo>
#include <mpi.h>

using namespace std;
using namespace arma;


// calculate and return the max eigenvalue of A' * A
double powIter (mat A);


//---------------------------------------------------------------------

int main () {

	// Initialize the MPI environment
	MPI_Init(NULL, NULL);
	int nProcs, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	srand(rank*12345);	// seed each proc

	
	

	// UNIT TESTS
	mat A = randu<mat>(5,5);
	cout << A.n_rows << "x" << A.n_cols << " matrix A: " << endl;
	cout << A << endl;
	cout << "Max singular value: \n" << max(svd(A)) << endl;
	cout << (powIter(A) - max(svd(A))) << endl;	
	return 0;
}

//---------------------------------------------------------------------

/*	
 *	Power Method : Parallel Implementation
 *	INPUT 
 *		- Matrix A
 *		-
 *	OUTPUT
 *		- The first sigma value of A
 *		- 
 *
 * 	METHOD
 * 		- Use the convergence of A' * A to find sigma^2
 * 		- 	
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

		// Use ALLREDUCE to sum each subsection of q together
		int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)

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


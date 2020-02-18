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
// 2/18/20
// MPI Incorporated into the code, but with runtime errors
// 	- ALL_REDUCE included
// 	- Additional variables added
// 	- MPI initialization and finalization added
//	- various fixme problems added
//
//
//

#include <iostream>
#include <armadillo>
#include <mpi.h>

using namespace std;
using namespace arma;


// calculate and return the max eigenvalue of A' * A
// implemented in parallel
double powIter (mat A, int &nProcs, int &rank);


//---------------------------------------------------------------------

int main () {

	// Initialize the MPI environment
	MPI_Init(NULL, NULL);
	int nProcs, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	srand(rank*12345);	// seed each proc
	

	// UNIT TESTS
	if (rank == 1) {
		mat A = randu<mat>(5,5);
		cout << A.n_rows << "x" << A.n_cols << " matrix A: " << endl;
		cout << A << endl;
		cout << "Max singular value: \n" << max(svd(A)) << endl;
		cout << (powIter(A, nProcs, rank) - max(svd(A))) << endl;	
	}
	
	MPI_Finalize();
	return 0;
}

//---------------------------------------------------------------------

/*	
 *	Power Method : Parallel Implementation
 *	INPUT 
 *		- Matrix A
 *		- nprocs and rank for PUs
 *	OUTPUT
 *		- The first sigma value of A
 * 	METHOD
 * 		- Use the convergence of A' * A to find sigma
 */
double powIter (mat A, int &nProcs, int &rank) { // FIXME Pass MAT A by reference?
	
	// FIXME NOTE: Armadillo stored column major format (potential time improvement?)
	
	int M = size(A,1);
	int N = size(A,2);

	// generate random tmpQ vector for each processor's submatrix operations
	// and a finalQ denoted for the original A matrix
	vec tmpQ(N);
	vec finalQ = randu<vec>(N);
	finalQ = finalQ / norm(finalQ, 2);
	cout << "TEST!!!!!!!!!!!!!!!" << endl;

	// prepare variables
	vec z;
	auto sigma = norm(finalQ,2); // FIXME norm of finalQ or tmpQ?
	auto s2 = sigma;
	double epsilon = 1.0;

	cout << "\nConvergence: " << endl;

	// converge to first sigma value of A
	int iter = 0;
	while (iter < 10 && epsilon > 0.00001) {
		z = A * finalQ;		// FIXME also allreduce for z?	
		cout << "z1: " << rank << "\n" << z << endl;
		z = z / norm(z,2);	// normalize z to manage large numbers
		cout << "z2: " << rank << "\n" << z << endl;
		tmpQ = A.t() * z;		

		// Use ALLREDUCE to sum each tmpQ together
		// Armadillo is high quality and can be passed through MPI - do not need arrays
		MPI_Allreduce(tmpQ.begin(), finalQ.begin(), N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		cout << "AFTER" << endl;

		sigma = norm(finalQ,2);	// set sigma here to prevent repetitive calculations
		finalQ = finalQ / sigma;
			
		epsilon = abs(sigma - s2)/(sigma);
		s2 = sigma;
		iter++;

		// debugging output
		cout << "Sigma:  " << sigma << endl;
	}

	return sigma;
}


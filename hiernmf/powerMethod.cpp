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
// 2/19/20
// MPI Fully functional
//

#include <iostream>
#include <armadillo>
#include <mpi.h>

using namespace std;
using namespace arma;


// calculate and return the max eigenvalue of A' * A
// implemented in parallel
double powIter (mat &localA, int &nProcs, int &rank);


//---------------------------------------------------------------------

int main () {

	// Initialize the MPI environment
	MPI_Init(NULL, NULL);
	int nProcs, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	srand(1);
	
	// Create A and partition for each PU's local matrices
	mat A = randu<mat>(4,4);
	if (rank == 0) {
		cout << A.n_rows << "x" << A.n_cols << " matrix A: " << endl;
		cout << A << endl;
		cout << "Max singular value: \n" << max(svd(A)) << endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	// Partition rows to PUs
	int M,N;
	M = size(A,0);
	N = size(A,1);
	int Pr = M / nProcs;
	int Pc = N / nProcs;
	int startRow = Pr * rank;
	int endRow = startRow + Pr;
	mat At = A.t();
	mat localA = At.cols(startRow, endRow-1);
	localA = localA.t();

	
	// power iteration:
	double sigmaSquared = powIter(localA,nProcs,rank);
	double sigma = max(svd(A));
	cout << abs(sigmaSquared - sigma*sigma) << endl;	


	MPI_Finalize();
	return 0;
}


/*	
 *	Power Method : Parallel Implementation
 *	INPUT 
 *		- Matrix A
 *		- nprocs and rank for PUs
 *	OUTPUT
 *		- The first sigma value of A squared
 */
double powIter (mat &localA, int &nProcs, int &rank) {
	
	int M = size(localA,0);
	int N = size(localA,1);

	// generate random localQ vector for each processor's submatrix operations
	// and a globalQ denoted for the original A matrix
	vec localQ(N);
	vec globalQ = randu<vec>(N);
	globalQ = globalQ / norm(globalQ, 2);

	// prepare variables
	vec z;
	auto sigma = norm(globalQ,2);
	auto s2 = sigma;
	double epsilon = 1.0;

	// DEBUG
	// cout << "\nConvergence: " << endl;

	// converge to first sigma value of A
	int iter = 0;
	while (iter < 10 && epsilon > 0.00001) {
		z = localA * globalQ;
		localQ = localA.t() * z;		

		// sum localQ into globalQ
		MPI_Allreduce(localQ.begin(), globalQ.begin(), N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		sigma = norm(globalQ,2);
		globalQ = globalQ / sigma;
			
		epsilon = abs(sigma - s2)/(sigma);
		s2 = sigma;
		iter++;

		// DEBUG
		// cout << "Sigma:  " << sigma << endl;
	}

	return sigma;
}


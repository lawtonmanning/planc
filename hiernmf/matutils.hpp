#include "common/utils.hpp"

namespace planc {
  template <class INPUTMATTYPE>
    void print(INPUTMATTYPE M, const char * name) {
      int size, rank;
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);

      MPI_Barrier(MPI_COMM_WORLD);
      if (rank == 0) {
        M.print(name);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      for (int i = 1; i < size; i++) {
        if (i == rank) {
          M.print();
        }
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }
  
  template <class INPUTMATTYPE>
  auto powIter (INPUTMATTYPE &A, int max_iter, double tol) -> double {

    int M = size(A,0);
    int N = size(A,1);

    // generate random localQ VECtor for each processor's submatrix operations
    // and a globalQ denoted for the original A matrix
    VEC localQ(N);
    VEC globalQ = arma::randu<VEC>(N);
    globalQ = globalQ / norm(globalQ, 2);

    // prepare variables
    VEC z;
    auto sigma = norm(globalQ,2);
    auto s2 = sigma;
    double epsilon = tol+1;

    // converge to first sigma value of A
    int iter = 0;
    while (iter < max_iter && epsilon > tol) {
      z = A * globalQ;
      localQ = z * A;		

      // sum localQ into globalQ
      MPI_Allreduce(localQ.begin(), globalQ.begin(), N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      sigma = norm(globalQ,2);
      globalQ = globalQ / sigma;

      epsilon = abs(sigma - s2)/(sigma);
      s2 = sigma;
      iter++;
    }

    return sigma;
  }

  VEC maxk(VEC X, int k) {
    VEC Xs = arma::sort(X, "descend");
    if (X.n_elem <= k) {
      return Xs;
    }
    return Xs.head(k);
  }

  UVEC maxk_idx(VEC X, int k) {
      UVEC Xi = arma::sort_index(X, "descend");
      if (X.n_elem <= k) {
        return Xi;
      }
      return Xi.head(k);
  }
}
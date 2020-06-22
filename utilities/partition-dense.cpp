#include <cstdio>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <string>
#include <iostream>
#include <cinttypes>
#include <armadillo>

/**
 * The dimension a particular rank holds out of 
 * the global dimension n across p processes.
 * @param[in] n is the global size of the tensor on that dimension
 * @param[in] p is the number of splits of n. 
 *            Typically number of processes on a particular mode
 * @param[in] r is the rank of the mpi process in that mode. fiber rank
 */
inline uint64_t itersplit(uint64_t n, int p, int r) {
  int split = (r < n % p) ? n / p + 1 : n / p;
  return split;
}

/**
 * Returns the start idx of the current rank r
 * for a global dimension n across p processes. 
 * @param[in] n is the size of the tensor on that dimension
 * @param[in] p is the number of splits of n.
 * @param[in] r is the rank of the mpi process in that mode. fiber rank
 */
inline uint64_t startidx(uint64_t n, int p, int r) {
  int rem = n % p;
  int idx =
      (r < rem) ? r * (n / p + 1) : (rem * (n / p + 1) + ((r - rem) * (n / p)));
  return idx;
}

/**
 * Returns the rank of the processor group
 * that owns the value at entry idx i for a
 * global dimension n across p processes.
 * @param[in] n is the size of the tensor on that dimension
 * @param[in] p is the number of splits of n
 * @param[in] i is the idx of the entry along that dimension
 */
inline int idxproc(int n, int p, int i) {
 int rem = n % p;
 int idx = 
     (i < rem * (n / p + 1)) ? i / (n / p + 1) : rem + (i - rem * (n / p + 1))/(n / p);
 return idx;
}



int main(int argc, char **argv) {
  printf("partition-dense began.\n");
  if (argc < 4) {
    printf("Usage: partition-dense [matrix-file-name] [row-proc-count] [col-proc-count]\n");
    return 0;
  }

  arma::mat X;
  X.load(argv[1]);
  int pr = atoi(argv[2]);
  int pc = atoi(argv[3]);
  uint64_t m = X.n_rows;
  uint64_t n = X.n_cols;

  printf("m::%" PRIu64 "::n::%" PRIu64 "::pr::%d::pc::%d\n",m,n,pr,pc);

  int idx = 0;
  for (int i = 0; i < pr; i++) {
    for (int j = 0; j < pc; j++) {
      printf("Writing the matrix for part %d...\n", idx);
      uint64_t m_start = startidx(m,pr,i);
      uint64_t n_start = startidx(n,pc,j);
      uint64_t m_end = m_start + itersplit(m,pr,i) - 1;
      uint64_t n_end = n_start + itersplit(n,pc,j) - 1;
      arma::mat x = X.submat(m_start,n_start , m_end,n_end);
      std::string out(argv[1]);
      out += std::to_string(idx++);
      x.save(out.c_str(),arma::arma_ascii);
    }
  }
}
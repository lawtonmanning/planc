#include <cstdio>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <string>
#include <iostream>
#include <cinttypes>

/**
 * The dimension a particular rank holds out of 
 * the global dimension n across p processes.
 * @param[in] n is the global size of the tensor on that dimension
 * @param[in] p is the number of splits of n. 
 *            Typically number of processes on a particular mode
 * @param[in] r is the rank of the mpi process in that mode. fiber rank
 */
inline int itersplit(int n, int p, int r) {
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
inline int startidx(int n, int p, int r) {
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
  printf("partition-uniform began.\n");
  if (argc < 4) {
    printf("Usage: partition-uniform [matrix-file-name] [row-proc-count] [col-proc-count]\n");
    return 0;
  }
  int rowProcCount = atoi(argv[2]);
  int colProcCount = atoi(argv[3]);
  int procCount = rowProcCount * colProcCount;
  printf("Creating a %d x %d uniform partition.\n", rowProcCount, colProcCount);

  printf("Opening input and output files...\n");
  FILE *file = fopen(argv[1], "r");
  if (file == NULL) {
    printf("Unable to open file %s.\n", argv[1]);
    return 0;
  }

  std::vector< std::vector<uint64_t> > procRowIdxs(procCount);
  std::vector< std::vector<uint64_t> > procColIdxs(procCount);
  std::vector< std::vector<double> > procVals(procCount);
  printf("Partitioning nonzeros and outputting to files...\n");
  int order;
  uint64_t nnz;
  uint64_t rowCount, colCount;
  fscanf(file, "%" SCNu64, &rowCount);
  fscanf(file, "%" SCNu64, &colCount);
  fscanf(file, "%" SCNu64, &nnz);
  printf("order::%d::nnz::%ju::rowcount::%ju::colcount::%ju\n", order, nnz, rowCount, colCount);
  for (uint64_t i = 0; i < nnz; i++) {
    uint64_t rowIdx, colIdx;
    double val;
    fscanf(file, "%" SCNu64, &rowIdx);
    fscanf(file, "%" SCNu64, &colIdx);
    fscanf(file, "%lf", &val);
    rowIdx--; colIdx--;
    uint64_t procRowIdx = idxproc(rowCount,rowProcCount,rowIdx);
    uint64_t procColIdx = idxproc(colCount,colProcCount,colIdx);
    if (procRowIdx >= rowProcCount || procColIdx >= colProcCount) { continue; }
    uint64_t procIdx = procRowIdx * colProcCount + procColIdx;
    uint64_t localRowIdx = rowIdx - startidx(rowCount,rowProcCount,procRowIdx);
    uint64_t localColIdx = colIdx - startidx(colCount,colProcCount,procColIdx);
    procRowIdxs[procIdx].push_back(localRowIdx);
    procColIdxs[procIdx].push_back(localColIdx);
    procVals[procIdx].push_back(val);
  }
  /*
  uint64_t rowsPerProc = rowCount / rowProcCount;
  uint64_t colsPerProc = colCount / colProcCount;
  printf("Assigning %ju rows and %ju columns per part.\n", rowsPerProc, colsPerProc);
  for (uint64_t i = 0; i < nnz; i++) {
    //if (i % 100000 == 0 && i > 0) { printf("Processing %ju.th nonzero...\n", i); }
    uint64_t rowIdx, colIdx;
    double val;
    fscanf(file, "%" SCNu64, &rowIdx);
    fscanf(file, "%" SCNu64, &colIdx);
    fscanf(file, "%lf", &val);
    rowIdx--; colIdx--;
    uint64_t procRowIdx = rowIdx / rowsPerProc;
    uint64_t procColIdx = colIdx / colsPerProc;
    if (procRowIdx >= rowProcCount || procColIdx >= colProcCount) { continue; } // Prune matrix
    uint64_t procIdx = procRowIdx * colProcCount + procColIdx;
    uint64_t localRowIdx = rowIdx % rowsPerProc;
    uint64_t localColIdx = colIdx % colsPerProc;
    procRowIdxs[procIdx].push_back(localRowIdx);
    procColIdxs[procIdx].push_back(localColIdx);
    procVals[procIdx].push_back(val);
  }
  */
  fclose(file);

  for (int i = 0; i < procCount; i++) {
    printf("Writing the matrix for part %d...\n", i);
    std::string outFileName(argv[1]);
    outFileName += std::to_string(i);
    file = fopen(outFileName.c_str(), "w");
    if (file == NULL) {
      printf("Unable to open file %s.\n", outFileName.c_str());
      return 0;
    }
    uint64_t m = itersplit(rowCount,rowProcCount,i/colProcCount);
    uint64_t n = itersplit(colCount,colProcCount,i%colProcCount);
    //printf("%" PRIu64 " %" PRIu64 " %" PRIu64 "\n", m, n, procVals[i].size());
    //fprintf(file, "%" PRIu64 " %" PRIu64 " %" PRIu64 "\n", m, n, procVals[i].size());
    auto &curRowIdxs = procRowIdxs[i];
    auto &curColIdxs = procColIdxs[i];
    auto &curVals = procVals[i];
    for (uint64_t j = 0; j < curRowIdxs.size(); j++) {
      fprintf(file, "%" PRIu64 " %" PRIu64 " %e\n", curRowIdxs[j], curColIdxs[j], curVals[j]);
    }
    fclose(file);
  }

  printf("partition-uniform finished.\n");
  return 0;
}


#include <cstdio>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <string>
#include <iostream>
#include <cinttypes>
#include <cmath>

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
inline int idxproc(uint64_t n, int p, uint64_t i) {
 int rem = n % p;
 int idx = 
     (i < rem * (n / p + 1)) ? i / (n / p + 1) : rem + (i - rem * (n / p + 1))/(n / p);
 return idx;
}

uint64_t range(uint64_t x[], int n) {
  uint64_t max = x[0];
  uint64_t min = x[0];
  for (int i = 1; i < n; i++) {
    if (x[i] > max) {
      max = x[i];
    }
    if (x[i] < min) {
      min = x[i];
    }
  }
  return max - min;
}

void sum(uint64_t sums[], std::vector< std::vector<uint64_t> > & counts, int p) {
  for (int i = 0; i < p; i++) {
    sums[i] = 0;
    for (uint64_t j = 0; j < counts[i].size(); j++) {
      sums[i] += counts[i][j];
    }
  }
}

inline uint64_t absdiff(uint64_t first, uint64_t second) {
  return (first < second) ? second - first : first - second;
}

void swap(std::vector< std::vector<uint64_t> > & counts, std::vector< std::vector<uint64_t> > & orders, uint64_t diff, int first, int second) {
  uint64_t imin = absdiff(diff, counts[first][0]);
  uint64_t iidx = 0;
  for (uint64_t i = 1; i < counts[first].size(); i++) {
    uint64_t new_min = absdiff(diff, counts[first][i]);
    if (new_min < imin) {
      imin = new_min;
      iidx = i;
    }
  }
  uint64_t x = counts[first][iidx];

  uint64_t jmin = absdiff(x / 2, counts[second][0]);
  uint64_t jidx = 0;
  for (uint64_t j = 1; j < counts[second].size(); j++) {
    uint64_t new_min = absdiff(diff, counts[second][j]);
    if (new_min < jmin) {
      jmin = new_min;
      jidx = j;
    }
  }
  uint64_t y = counts[second][jidx];

  counts[first][iidx] = y;
  counts[second][jidx] = x;

  uint64_t tmp = orders[first][iidx];
  orders[first][iidx] = orders[second][jidx];
  orders[second][jidx] = tmp;
}

void entropy(uint64_t counts[], uint64_t orders[], uint64_t n, int p) {
  std::vector< std::vector<uint64_t> > pcounts(p);
  std::vector< std::vector<uint64_t> > porders(p);
  uint64_t sums[p];
  for (int i = 0; i < p; i++) {
    uint64_t start = startidx(n,p,i);
    uint64_t end = start + itersplit(n,p,i);
    for (uint64_t j = start; j < end; j++) {
      pcounts[i].push_back(counts[j]);
      porders[i].push_back(orders[j]);
    }
  }
  sum(sums,pcounts,p);
  for (int i = 0; i < p; i++) {
    printf("proc::%d::nnz::%ju\n",i,sums[i]);
  }
  uint64_t old_range = range(sums,p);
  uint64_t new_range = old_range-1;
  std::vector< std::vector<uint64_t> > old_pcounts = pcounts;
  std::vector< std::vector<uint64_t> > old_porders = porders;
  while (new_range < old_range) {
    old_pcounts = pcounts;
    old_porders = porders;
    for (int i = 0; i < p; i++) {
      for (int j = i+1; j < p; j++) {
        sum(sums,pcounts,p);
        uint64_t diff = (sums[i] > sums[j]) ? sums[i] - sums[j] : sums[j] - sums[i];
        if (diff == 0) { continue; }
        else if (sums[i] > sums[j]) {
          swap(pcounts,porders,diff,i,j);
        }
        else {
          swap(pcounts,porders,diff,j,i);
        }
      }
    }
    old_range = new_range;
    sum(sums,pcounts,p);
    new_range = range(sums,p);
    printf("%ju\n",new_range);
  }
  pcounts = old_pcounts;
  porders = old_porders;

  sum(sums,pcounts,p);
  new_range = range(sums,p);
  printf("%ju\n",new_range);
  for (int i = 0; i < p; i++) {
    printf("proc::%d::nnz::%ju\n",i,sums[i]);
  }

  for (int i = 0; i < p; i++) {
    uint64_t start = startidx(n,p,i);
    uint64_t end = start + itersplit(n,p,i);
    uint64_t idx = 0;
    for (uint64_t j = start; j < end; j++) {
      orders[j] = porders[i][idx];
      idx++;
    }
  }

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
    printf("Unable to open file %s\n", argv[1]);
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
  
  uint64_t rowCounts[rowCount];
  uint64_t colCounts[colCount];
  uint64_t rowOrders[rowCount];
  uint64_t colOrders[colCount];
  
  for (uint64_t i = 0; i < rowCount; i++) {
    rowCounts[i] = 0;
    rowOrders[i] = i;
  }
  for (uint64_t i = 0; i < colCount; i++) {
    colCounts[i] = 0;
    colOrders[i] = i;
  }
  
  std::vector<uint64_t> rowIdxs(nnz);
  std::vector<uint64_t> colIdxs(nnz);
  std::vector<double> vals(nnz);
  for (uint64_t i = 0; i < nnz; i++) {
    uint64_t rowIdx, colIdx;
    double val;
    fscanf(file, "%" SCNu64, &rowIdx);
    fscanf(file, "%" SCNu64, &colIdx);
    fscanf(file, "%lf", &val);
    rowIdx--; colIdx--;
    rowCounts[rowIdx]++;
    colCounts[colIdx]++;
    rowIdxs[i] = rowIdx;
    colIdxs[i] = colIdx;
    vals[i] = val;
  }

  entropy(rowCounts, rowOrders, rowCount, rowProcCount);

  for (uint64_t i = 0; i < nnz; i++) {
    uint64_t rowIdx, colIdx;
    double val;
    rowIdx = rowOrders[rowIdxs[i]];
    colIdx = colOrders[colIdxs[i]];
    val = vals[i];
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


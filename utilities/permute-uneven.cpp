#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <ctime>
#include <string>
#include <iostream>
#include <cinttypes>

int main(int argc, char **argv) {
  printf("permute-uneven began.\n");
  if (argc != 5) {
    printf("Usage: permute-uneven [matrix-file-name] [output-file-name] [row-idxs-name] [col-idxs-name]\n");
    return 0;
  }

  printf("Opening input and output files...\n");
  FILE * ifile = fopen(argv[1], "r");
  if (ifile == NULL) {
    printf("Unable to open file %s.\n", argv[1]);
    return 0;
  }
  FILE * ofile = fopen(argv[2], "w");
  if (ofile == NULL) {
    printf("Unable to open file %s.\n", argv[2]);
    return 0;
  }

  uint64_t nnz;
  uint64_t rowCount, colCount;
  fscanf(ifile, "%" SCNu64, &rowCount);
  fscanf(ifile, "%" SCNu64, &colCount);
  fscanf(ifile, "%" SCNu64, &nnz);
  printf("nnz::%ju::rowcount::%ju::colcount::%ju\n", nnz, rowCount, colCount);
  fprintf(ofile, "%" PRIu64 " %" PRIu64 " %" PRIu64 "\n", rowCount, colCount, nnz);
  std::vector<uint64_t> rowIdxs(rowCount);
  std::vector<uint64_t> colIdxs(colCount);
  for (uint64_t i = 0; i < rowCount; i++) {
    rowIdxs[i] = i;
  }
  for (uint64_t i = 0; i < colCount; i++) {
    colIdxs[i] = i;
  }
  std::random_shuffle(rowIdxs.begin(), rowIdxs.end());
  std::random_shuffle(colIdxs.begin(), colIdxs.end());
  for (uint64_t i = 0; i < nnz; i++) {
    uint64_t rowIdx, colIdx;
    double val;
    fscanf(ifile, "%" SCNu64, &rowIdx);
    fscanf(ifile, "%" SCNu64, &colIdx);
    fscanf(ifile, "%lf", &val);
    rowIdx--; colIdx--;

    fprintf(ofile, "%" PRIu64 " %" PRIu64 " %e\n", rowIdxs[rowIdx]+1, colIdxs[colIdx]+1, val);
  }
  fclose(ifile);
  fclose(ofile);

  FILE * rfile = fopen(argv[3], "w");
  if (rfile == NULL) {
    printf("Unable to open file %s.\n", argv[3]);
    return 0;
  }
  for (uint64_t i = 0; i < rowCount; i++) {
    if (i == 0) {
      fprintf(rfile, "%" PRIu64, rowIdxs[i]);
    }
    else {
      fprintf(rfile, "\n%" PRIu64, rowIdxs[i]);
    }
  }
  fclose(rfile);
  FILE * cfile = fopen(argv[4], "w");
  if (cfile == NULL) {
    printf("Unable to open file %s.\n", argv[4]);
    return 0;
  }
  for (uint64_t i = 0; i < colCount; i++) {
    if (i == 0) {
      fprintf(cfile, "%" PRIu64 , colIdxs[i]);
    }
    else {
      fprintf(cfile, "\n%" PRIu64 , colIdxs[i]);
    }
  }
  fclose(cfile);


  printf("permute-uneven finished.\n");
  return 0;
}


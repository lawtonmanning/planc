
using namespace planc;

template <class INPUTMATTPYE>
void print(INPUTMATTYPE M, const char * name) {
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
  fprintf(name, "%s:", name);
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


  

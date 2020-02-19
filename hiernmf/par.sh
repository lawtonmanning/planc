# compiling armadillo for c++ requires an additional parameter
# mpi compiled using mpicxx
mpicxx powerMethod.cpp -o parPowerMethod -larmadillo

mpirun -n 4 ./parPowerMethod

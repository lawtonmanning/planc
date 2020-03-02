#! /bin/bash

cd ../build/dense_distnmf
make || exit

echo "Even Division (32x32 on 4x4) > even"
mpirun -np 16 ./dense_distnmf -a 7 -i rand_lowrank -d "32 32" -p "4 4" -k 2 -t 2000 -e 1 > ~/planc/distnmf/even
grep "relerr" ~/planc/distnmf/even | tail -1

echo "Odd Division (33x33 on 4x4) > odd"
mpirun -np 16 ./dense_distnmf -a 7 -i rand_lowrank -d "33 33" -p "4 4" -k 2 -t 2000 -e 1 > ~/planc/distnmf/odd
grep "relerr" ~/planc/distnmf/odd | tail -1

echo "Odd 1DRow (33x33 on 1x16) > 1DRow"
mpirun -np 16 ./dense_distnmf -a 7 -i rand_lowrank -d "33 33" -p "1 16" -k 2 -t 2000 -e 1 > ~/planc/distnmf/1DRow
grep "relerr" ~/planc/distnmf/1DRow | tail -1

echo "Odd 1DCol (33x33 on 16x1) > 1DCol"
mpirun -np 16 ./dense_distnmf -a 7 -i rand_lowrank -d "33 33" -p "16 1" -k 2 -t 2000 -e 1 > ~/planc/distnmf/1DCol
grep "relerr" ~/planc/distnmf/1DCol | tail -1

echo "Odd C<P (33x12 on 2x8) > smallC"
mpirun -np 16 ./dense_distnmf -a 7 -i rand_lowrank -d "33 12" -p "2 8" -k 2 -t 2000 -e 1 > ~/planc/distnmf/smallC
grep "relerr" ~/planc/distnmf/smallC | tail -1

echo "Odd R<P (12x33 on 8x2) > smallR"
mpirun -np 16 ./dense_distnmf -a 7 -i rand_lowrank -d "12 33" -p "8 2" -k 2 -t 2000 -e 1 > ~/planc/distnmf/smallR
grep "relerr" ~/planc/distnmf/smallR | tail -1

echo "Odd Large (1231x1423 on 2x8) > large"
mpirun -np 16 ./dense_distnmf -a 7 -i rand_lowrank -d "1231 1423" -p "2 8" -k 2 -t 2000 -e 1 > ~/planc/distnmf/large
grep "relerr" ~/planc/distnmf/large | tail -1

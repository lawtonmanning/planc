#!/bin/bash
#SBATCH -A CSC345
#SBATCH -J R2NMF
#SBATCH -N 1
#SBATCH -t 1:00:00
#SBATCH -p gpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mannlg15@wfu.edu
#SBATCH --get-user-env

module list
cd ~/planc/build-rhea
cmake ..
make || exit

pi=($(seq 1 56))
pj=($(seq 1 56))

for i in "${pi[@]}"
do
  for j in "${pj[@]}"
  do
    p=$(($i * $j))
    if [ "$p" -le 56 ]; then
      srun -n $p ./dense_distnmf/dense_distnmf -a 3 -i rand_lowrank -d "4096 4096" -p "$i $j" -k 2 -t 2000 -e 0 >> ~/planc/build-rhea/proc$p
    fi
  done
done

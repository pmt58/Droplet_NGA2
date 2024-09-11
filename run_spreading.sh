#!/bin/bash
#CA=(10 50)
#($(seq 10 50))
#OH=(1e-3 3e-3 1e-2 3e-2 1e-1 3e-1 1e+0 3e+0 1e+1 3e+1 1e+2)
Beta=(0.1 0.2 0.4 0.6 0.8 1)
#for j in ${!CA[@]}; do
#   for i in ${!OH[@]}; do
#      echo "Running a simulation with Oh=${OH[$i]} and CA=${CA[$j]}"
#      mpiexec -n 4 ./nga.dp.gnu.opt.mpi.exe -i input  --Oh=${OH[$i]} --CA=${CA[$j]} >> result_${CA[$j]}.txt
#   done
#done

# Just one test run with Oh=1e-3
for j in ${!Beta[@]}; do
   echo "Running a simulation with Beta=${Beta[$j]}"
   mpiexec -n 4 ./nga.dp.gnu.opt.mpi.exe -i input --Beta=${Beta[$j]} >> result.txt
   mkdir ./monitor_${Beta[$j]}
   mv ./monitor/* ./monitor_${Beta[$j]}/
done

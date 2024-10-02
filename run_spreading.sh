#!/bin/bash

#CA=(10 50)
#($(seq 10 50))
OH=(1e-3 1e-1)
CL=(0 1)
#Beta=(0.06)
#G=(0)

# Iteratively test: High/Low Oh : CL model type : Mesh Resolution
for o in ${!OH[@]}; do
   for c in ${!CL[@]}; do
         echo "Running a simulation with Oh=${OH[$o]} CL=${CL[$c]}"
         mpiexec -n 4 ./nga.dp.gnu.opt.mpi.exe -i input --Oh=${OH[$o]} --CLsolver=${CL[$c]} >> result.txt
         mkdir ./data_${o}_${CL[$c]}
         mv ./monitor/* ./data_${o}_${CL[$c]}/
   done
done

#!/bin/bash

#CA=(10 50)
#($(seq 10 50))
OH=(0.3) #(1e-3 1e-1)
CL=(1) #1 2)
Nele=(20 30 50 70)
#Beta=(0.06)
#G=(0)

# Iteratively test: High/Low Oh : CL model type : Mesh Resolution
for o in ${!OH[@]}; do
   for c in ${!CL[@]}; do
      for n in ${!Nele[@]}; do
         echo "Running a simulation with Oh=${OH[$o]} CL=${CL[$c]} N=${Nele[$n]}"
         mpiexec -n 4 ./nga.dp.gnu.opt.mpi.exe -i input --Oh=${OH[$o]} --CLsolver=${CL[$c]} --nx=${Nele[$n]} --ny=${Nele[$n]} --nz=${Nele[$n]} >> result.txt
         mkdir ./data1_${CL[$c]}_${o}_${n}
         mv ./monitor/* ./data1_${CL[$c]}_${o}_${n}/
      done
   done
done

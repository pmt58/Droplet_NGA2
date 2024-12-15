#!/bin/bash

#CA=(10 50)
#($(seq 10 50))
OH=(0.05) #(1e-3 1e-1)
CL=(1)
Nele=(32)
Lamb=(0.1 1 10)
#G=(0)

# Iteratively test: High/Low Oh : CL model type : Mesh Resolution
for o in ${!OH[@]}; do
   for c in ${!CL[@]}; do
      for n in ${!Nele[@]}; do
         for l in ${!Lamb[@]}; do
            echo "Running a simulation with Oh=${OH[$o]} CL=${CL[$c]} N=${Nele[$n]}"
            mpiexec -n 4 ./nga.dp.gnu.opt.mpi.exe -i input --Lambda=${Lamb[$l]} --Oh=${OH[$o]} --CLsolver=${CL[$c]} --nx=${Nele[$n]} --ny=${Nele[$n]} --nz=${Nele[$n]} >> result.txt
            mkdir ./data4_${CL[$c]}_${l}_${n}
            mv ./monitor/* ./data4_${CL[$c]}_${l}_${n}/
         done
      done
   done
done

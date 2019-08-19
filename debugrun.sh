#!/bin/bash 

compiler="gfortran"
options="-fbounds-check -Wall -Wno-uninitialized -Wno-unused-function -Wno-unused-dummy-argument -fimplicit-none -fcheck=all -g -pedantic"
source="ModIntVec.f90 ModVec.f90 ModMtrx.f90 ModSparse.f90 ModAppr.f90 ModRecon.f90 Main.f90"
object="ModIntVec.o ModVec.o ModMtrx.o ModSparse.o ModAppr.o ModRecon.o Main.o"
program="Main_gnu.exe"

mkdir -p obj_gnu;
for file in $source; do
  $compiler -J./obj_gnu -I./obj_gnu $options -c src/$file
  if [[ $? != 0 ]]; then
    exit
  fi
done
$compiler -I./obj_gnu -o $program $object -llapack -lblas
mv $object obj_gnu
./$program

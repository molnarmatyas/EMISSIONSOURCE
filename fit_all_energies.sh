#!/bin/bash 

# Redirect both stdout and stderr to a file
exec > fitall.log 2>&1

# Array of inputs 
energies=("7.7GeV" "9.2GeV" "19.6GeV") # "27GeV" does not work yet 

for energy in "${!energies[@]}"
do
  echo "Running processing and fit for energy: ${energies[$energy]}"
  root.exe -b -q "eposSpatialDistKT.cpp(100, 100, 1, 0, \"${energies[$energy]}\")"
#  root.exe -b -q distPoints.cpp
  for i in {1..4}
  do
    echo "Running for KT bin ${i}"
    ./levy_fit.exe $i ${energies[$energy]}
  done
done


#echo "Running processing and fit for energy: 27GeV"
#root.exe -b -q 'eposSpatialDistKT.cpp(100, 50, 1, 0, "27GeV")'
#root.exe -b -q distPoints.cpp
#for i in {1..5}
#do
#  ./levy_fit.exe $i 27GeV
#done

echo "Running processing and fit for energy: 200GeV"
root.exe -b -q 'eposSpatialDistKT.cpp(100, 25, 1, 0, "200GeV")'
#root.exe -b -q distPoints.cpp
for i in {1..4}
do
  ./levy_fit.exe $i 200GeV
done

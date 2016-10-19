#!/bin/bash


# Read command line argument
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
-p|--pathToQt5)
PATHTOQT="$2"
shift # past argument
;;
*)
# unknown option
;;
esac
shift # past argument or value
done


if [${PATHTOQT} -eq ""]
then
  echo "Please specify the -p argument (path to the cmake folder of Qt5)"
else

cd src
# Cmake to configure the project
cmake CMakeLists.txt -DPATHTOQT=${PATHTOQT}
# Compile the source code
make
# Copy the programs into the bin folder
mkdir ../bin
cp FlatFieldCorrection/FlatFieldCorrection ../bin/FlatFieldCorrection
cp PreMosa/PreMosa ../bin/PreMosa
cp Projection/Projection ../bin/Projection
cp ProjectionUsingHeightMap/ProjectionUsingHeightMap ../bin/ProjectionUsingHeightMap

cd ..

echo -e "\n\n---> Successful building! Executables under ./bin/"

fi
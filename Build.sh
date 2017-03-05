#!/bin/bash

cd src
# Cmake to configure the project
cmake CMakeLists.txt
# Compile the source code
make
# Copy the programs into the bin folder
mkdir ../bin
cp FlatFieldCorrection/FlatFieldCorrection ../bin/FlatFieldCorrection
cp PreMosa/PreMosa ../bin/PreMosa
cp SurfaceExtraction/SurfaceExtraction ../bin/SurfaceExtraction
cp ExtendedSurfaceExtraction/ExtendedSurfaceExtraction ../bin/ExtendedSurfaceExtraction

cd ..

echo -e "\n\n---> Successful building! Executables under ./bin/"
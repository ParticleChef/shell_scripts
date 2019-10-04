#!/bin/bash
mkdir -p Result_files

mkdir -p Result_files/${1}
mkdir -p Result_files/${1}/txtfiles
mkdir -p Result_files/${1}/codes
mkdir -p Result_files/${1}/reports

cp MET* Result_files/${1}/codes/.
cp met/solution_test/syn/report/* Result_files/${1}/reports
cp met/solution_test/csim/build/*.txt Result_files/${1}/txtfiles




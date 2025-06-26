#! /bin/bash

# shell script to process SAXS data with USToo programs
# and generate a Guinier plot
# accepts a single input argument, the name of the saxs data file
# uses directed input and output with pipes
# suppresses intermediate plots 
# and saves only the final output file from saxsDeSmear

dataFileName=$1
bufferFileName="RNAseA_buffer.pdh"
rootName=${dataFileName%%.pdh}
ext="_diff_bin_dsm.pdh"
outFileName="$rootName$ext"
saxsSubtract $bufferFileName $dataFileName -t none --so | \
saxsBinData -n 5 --si -t none | \
saxsDeSmear desmearParam.txt --si > $outFileName
saxsGuinier 0.03 0.06 $outFileName

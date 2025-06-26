#! /bin/bash

# shell script to process SAXS data with USToo programs
# and generate a Guinier plot
# uses directed input and output with pipes
# suppresses intermediate plots 
# and saves only the final output file from saxsDeSmear

saxsSubtract RNAseA_buffer.pdh RNAseA.pdh -t none --so | \
saxsBinData -n 5 --si -t none | \
saxsDeSmear desmearParam.txt --si > RNAseA_diff_bin_dsm.txt 
saxsGuinier 0.03 0.06 RNAseA_diff_bin_dsm.txt

#! /bin/bash

# shell script to process SAXS data with USToo programs
# and generate a Guinier plot
# shows plots at each step and saves intermediate files
saxsSubtract RNAseA_buffer.pdh RNAseA.pdh
saxsBinData -n 5 RNAseA_diff.pdh
saxsDeSmear desmearParam.txt RNAseA_diff_bin.pdh
saxsGuinier 0.03 0.06 RNAseA_diff_bin_dsm.pdh

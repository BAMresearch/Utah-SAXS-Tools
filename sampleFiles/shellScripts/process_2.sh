#! /bin/bash

# shell script to process SAXS data with USToo programs
# and generate a Guinier plot
# supresses intermediate plots, but saves intermediate files

saxsSubtract RNAseA_buffer.pdh RNAseA.pdh -t none
saxsBinData -n 5 RNAseA_diff.pdh -t none
saxsDeSmear desmearParam.txt RNAseA_diff_bin.pdh
saxsGuinier 0.03 0.06 RNAseA_diff_bin_dsm.pdh

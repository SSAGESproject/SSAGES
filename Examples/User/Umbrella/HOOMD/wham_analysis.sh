#!/bin/bash

# Usage: wham [P|Ppi|Pval] hist_min hist_max num_bins tol temperature numpad \
#        metadatafile freefile [num_MC_trials randSeed]
./wham Ppi -3.14159 3.14159 100 0.0001 300 1 wham_metafile wham_freefile > wham_output

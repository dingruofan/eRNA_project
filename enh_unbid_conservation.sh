#!/bin/bash
bwtool agg 2000:2000 ./enh_find_1/result.bed ../../phastCons46/phastCons46.bw  /dev/stdout > ./enh_find_1/result.mean_signal
now_path=`pwd`
Rscript $2/unbid_conservation_plot.R $now_path
feh $now_path/result.mean_signal.png


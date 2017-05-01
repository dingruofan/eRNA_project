cd /home/ding/all_cmd/script/enh_statistics
head ./enh_all_cluster1.bed -n 30
linkBed -db hg19 -i ./enh_all_cluster1.bed >./a.html

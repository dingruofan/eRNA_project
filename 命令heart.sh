#建立目录：
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
mkdir -p ./bed/H3K27ac ./bed/H3K4me1 ./bed/H3K4me3
mkdir -p ./bg
mkdir -p ./enh_find

#解压缩，可以手动。
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
gunzip -c GSM*H3K27ac*.bed.gz > H3K27ac.bed
gunzip -c GSM*H3K4me1*.bed.gz > H3K4me1.bed
gunzip -c GSM*H3K4me3*.bed.gz > H3K4me3.bed
gunzip -c GSM*Input*.bed.gz > chip_input.bed


#使用macs14获得组蛋白修饰位点:
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart/bed/H3K4me1
macs14 -t ../../H3K4me1.bed -c ../../chip_input.bed -n H3K4me1 -f BED -g hs --nomodel --shiftsize=73  -S -B 

cd /media/ding/000B49000006264C/eRNA_project/histone/Heart/bed/H3K4me3
macs14 -t ../../H3K4me3.bed -c ../../chip_input.bed -n H3K4me3 -f BED -g hs --nomodel --shiftsize=73  -S -B  

cd /media/ding/000B49000006264C/eRNA_project/histone/Heart/bed/H3K27ac
macs14 -t ../../H3K27ac.bed -c ../../chip_input.bed -n H3K27ac -f BED -g hs --nomodel --shiftsize=73  -S -B 

#删除bed文件来节约空间：
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
rm H3K27ac.bed H3K4me1.bed H3K4me3.bed chip_input.bed

#将bedgraph.gz转为bigwig文件：
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart/
/home/ding/tool/wigToBigWig ./bed/H3K4me3/H3K4me3_MACS_bedGraph/treat/H3K4me3_treat_afterfiting_all.bdg.gz /home/ding/tool/hg19.chrom.sizes  ./bg/H3K4me3.bw -clip
/home/ding/tool/wigToBigWig ./bed/H3K4me1/H3K4me1_MACS_bedGraph/treat/H3K4me1_treat_afterfiting_all.bdg.gz /home/ding/tool/hg19.chrom.sizes  ./bg/H3K4me1.bw -clip
/home/ding/tool/wigToBigWig ./bed/H3K27ac/H3K27ac_MACS_bedGraph/treat/H3K27ac_treat_afterfiting_all.bdg.gz /home/ding/tool/hg19.chrom.sizes  ./bg/H3K27ac.bw -clip


#查看是否排序：(都已经排序了)
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
bedextract --list-chr ./bed/H3K4me1/H3K4me1_peaks.bed
bedextract --list-chr ./bed/H3K27ac/H3K27ac_peaks.bed
bedextract --list-chr ./bed/H3K4me3/H3K4me3_peaks.bed



#######################################################################
#获取每一个组蛋白信号的总信号数目：
#######################################################################
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
bwtool summary ../../genome_wide.bed ./bg/H3K4me3.bw /dev/stdout -with-sum |awk 'BEGIN{a=0}(a=a+$NF){}END{print a}' - >./bg/H3K4me3_sigsum
bwtool summary ../../genome_wide.bed ./bg/H3K4me1.bw /dev/stdout -with-sum |awk 'BEGIN{a=0}(a=a+$NF){}END{print a}' - >./bg/H3K4me1_sigsum
bwtool summary ../../genome_wide.bed ./bg/H3K27ac.bw /dev/stdout -with-sum |awk 'BEGIN{a=0}(a=a+$NF){}END{print a}' - >./bg/H3K27ac_sigsum
bwtool summary ../../genome_wide.bed ./bg/mRNA.bw /dev/stdout -with-sum |awk 'BEGIN{a=0}(a=a+$NF){}END{print a}' - >./bg/mRNA_sigsum



#######################################################################
#启动子附近的组蛋白修饰：(太慢！)
#######################################################################
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
bwtool agg 3000:3000 ../../gencode/protein_lncRNA.ss ./bg/H3K4me1.bw,./bg/H3K4me3.bw,./bg/H3K27ac.bw -firstbase  /dev/stdout > ./enh_find/result.mean_signal



#######################################################################
#CAGE的enhancer预测的enhancer周围的组蛋白特征
#######################################################################

cd /media/ding/000B49000006264C/eRNA_project/fantom/Heart_CAGE
/home/ding/tool/enhancers-master/scripts/bidir_enhancers -d 0.8  -f ./ctss.path -o ./cage_enh  

#查看bed是否排序:(没有排序：)
cd /media/ding/000B49000006264C/eRNA_project/fantom/Heart_CAGE
bedextract --list-chr ./cage_enh/enhancers.bed
sort-bed ./cage_enh/enhancers.bed>./cage_enh/enhancers_sorted.bed
sort-bed ./cage_enh/bidir.pairs.bed>./cage_enh/bidir.pairs_sorted.bed


#获取双向转录的转录起始位点：
cd /media/ding/000B49000006264C/eRNA_project/fantom/Heart_CAGE
awk 'BEGIN{FS="\t"}split($11,a,","){print $1"\t"$2+a[1]"\t"$3-a[2]"\t"$4}' ./cage_enh/bidir.pairs_sorted.bed |bedops --not-element-of 1 - ../../gencode/protein_coding_gene_up1000 >./cage_enh/bidir.pairs_sorted_TSS.bed

#可以使用R做距离enhancer中心的密度图：

cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
closest-features --closest  --delim '\t' --dist ../../fantom/Heart_CAGE/cage_enh/enhancers_sorted.bed ./bed/H3K27ac/H3K27ac_peaks.bed |awk '($NF<3000&&$NF>-3000){print $NF}' -  >./enh_find/27ac.dist
closest-features --closest  --delim '\t' --dist ../../fantom/Heart_CAGE/cage_enh/enhancers_sorted.bed ./bed/H3K4me1/H3K4me1_peaks.bed |awk '($NF<3000&&$NF>-3000){print $NF}' -  >./enh_find/me1.dist
closest-features --closest  --delim '\t' --dist ../../fantom/Heart_CAGE/cage_enh/enhancers_sorted.bed ./bed/H3K4me3/H3K4me3_peaks.bed |awk '($NF<3000&&$NF>-3000){print $NF}' -  >./enh_find/me3.dist

#做组蛋白信号图：
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
bwtool agg 1000:1000 ../../fantom/Heart_CAGE/cage_enh/enhancers_sorted.bed ./bg/H3K4me1.bw,./bg/H3K4me3.bw,./bg/H3K27ac.bw /dev/stdout > ./enh_find/result.mean_signal



#######################################################################
#查看Dnase周围的组蛋白的特征：
#######################################################################
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
closest-features --closest  --delim '\t' --dist ../../DNase/Heart/Dnase_ENCODE.bed ./bed/H3K27ac/H3K27ac_peaks.bed |awk '($NF<3000&&$NF>-3000){print $NF}' -  >./enh_find/27ac.dist
closest-features --closest  --delim '\t' --dist ../../DNase/Heart/Dnase_ENCODE.bed ./bed/H3K4me1/H3K4me1_peaks.bed |awk '($NF<3000&&$NF>-3000){print $NF}' -  >./enh_find/me1.dist
closest-features --closest  --delim '\t' --dist ../../DNase/Heart/Dnase_ENCODE.bed ./bed/H3K4me3/H3K4me3_peaks.bed |awk '($NF<3000&&$NF>-3000){print $NF}' -  >./enh_find/me3.dist

cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
bwtool agg 3000:3000 ../../DNase/Heart/Dnase_ENCODE.bed ./bg/H3K4me1.bw,./bg/H3K4me3.bw,./bg/H3K27ac.bw -firstbase  /dev/stdout > ./enh_find/result.mean_signal



#######################################################################
#新思路2：筛选Dnase位点2000bp范围内组蛋白修饰的峰值summits，然而并没有summits
#######################################################################
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
closest-features --closest  --delim '\t' --dist ../../DNase/Heart/Dnase_ENCODE.bed ./bed/H3K27ac/H3K27ac_summits.bed |awk '($NF<2000&&$NF>-2000){print $1,$2,$3,$4,$NF,$6}' - |sed  's/ /\t/g' - >./enh_find/a
closest-features --closest  --delim '\t' --dist ./enh_find/a ./bed/H3K4me1/H3K4me1_summits.bed |awk '($NF<2000&&$NF>-2000){print $1,$2,$3,$4,$NF,$6}' - |sed  's/ /\t/g' - >./enh_find/b

#去掉蛋白编码区：
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
bedops --not-element-of 1 ./enh_find/b ../../gencode/protein_coding_gene_up1000 >./enh_find/c

#H3K4me3的处理：然而并没有处理：($NF>1||$NF<-1)
closest-features --closest  --delim '\t' --dist ./enh_find/c ./bed/H3K4me3/H3K4me3_summits.bed |awk '{print $1,$2,$3,$4,$NF,$6}' - |sed  's/ /\t/g' - >./enh_find/result.bed
wc -l ./enh_find/result.bed

#做蛋白信号图：
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
bwtool agg 3000:3000 ./enh_find/result.bed ./bg/H3K4me1.bw,./bg/H3K4me3.bw,./bg/H3K27ac.bw  /dev/stdout > ./enh_find/result.mean_signal

#######################################################################
#与CAGE定义的双向转录取交集：
#组蛋白定义的enhancer中心result.bed上下游扩增2kb:
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
awk '{print $1"\t"$2-2000"\t"$3+2000}' ./enh_find/result.bed  >./enh_find/result_2kb.bed

#2kb内有双向转录TSS的增强子：(bidir.pairs_sorted_TSS.bed)
#一个dnase 2kb范围内可能存在多个bidir，选择离他最近的那一个。
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
bedops --element-of ../../fantom/Heart_CAGE/cage_enh/bidir.pairs_sorted_TSS.bed ./enh_find/result_2kb.bed  |closest-features --closest  --delim '\t' --dist -  ./enh_find/result.bed >./enh_find/result_2kb_bidir_all.bed 
 
wc -l ./enh_find/result_2kb_bidir_all.bed 
#到此得到的./enh_find/result_2kb_bidir_all.bed 包含cage定义的双向转录为点（前4列，第四列是完整的），以及5,6,7列是与之最近的2000bp内的Dnase位点。

#处理一下列的位置并重命名：
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
awk '{print $5"\t"$6"\t"$7"\t"$4"\t"$11}'  ./enh_find/result_2kb_bidir_all.bed > ./enh_find/Dnase_histone_bidir_enh

#######################################################################

#获得2000bp内的转录因子及其结合位点。
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
closest-features --closest  --delim '\t' --dist ../../TFBS/TFBS_sorted.bed ./enh_find/Dnase_histone_bidir_enh |awk '($NF<=2000&&$NF>=-2000){print }' - |sed  's/ /\t/g' - >./enh_find/Dnase_histone_bidir_enh_TFBS
#列说明：前5列为转录因子结合为点信息,6,7,8列为Dnase结合为点，9列为bidir位点，10列为Dnase与bidir距离，11列为Dnase与TFBS距离。

#######################################################################
#添加TF表大量到最后一列：(表达量文件../../TF_exp/Heart/TF_exp)
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
awk 'BEGIN{FS="\t"}{if(NR==FNR){a[$2]=$3;next} if(NR>FNR){print $0"\t"a[$4]}}'  ../../TF_exp/Heart/TF_exp ./enh_find/Dnase_histone_bidir_enh_TFBS >./enh_find/Dnase_histone_bidir_enh_TFBS_exp 
wc -l ./enh_find/Dnase_histone_bidir_enh_TFBS_exp 
wc -l ./enh_find/Dnase_histone_bidir_enh

#######################################################################
#研究vista 增强子与我预测出来的
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
closest-features --closest  --delim '\t' --dist ../../vista_enhancer/heart/vista_enh.bed  ./enh_find/result.bed | awk '{print $NF}' ->./enh_find/vista_result.dist

cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
closest-features --closest  --delim '\t' --dist ../../vista_enhancer/heart/vista_enh.bed  ./enh_find/Dnase_histone_bidir_enh | awk '{print $NF}' ->./enh_find/vista_dnase_bidir.dist

#######################################################################

#可以使用R做距离enhancer中心的密度图：

cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
closest-features --closest  --delim '\t' --dist ./enh_find/result_cage.dist ./bed/H3K27ac/H3K27ac_peaks.bed |awk '($NF<2000&&$NF>-2000){print $NF}' -  >./enh_find/27ac.dist
closest-features --closest  --delim '\t' --dist ./enh_find/result_cage.dist ./bed/H3K4me1/H3K4me1_peaks.bed |awk '($NF<2000&&$NF>-2000){print $NF}' -  >./enh_find/me1.dist
closest-features --closest  --delim '\t' --dist ./enh_find/result_cage.dist ./bed/H3K4me3/H3K4me3_peaks.bed |awk '($NF<2000&&$NF>-2000){print $NF}' -  >./enh_find/me3.dist




##########################################################################
##转录起始位点附近的组蛋白和Dnase修饰：
#######################################################################
cd /media/ding/000B49000006264C/eRNA_project/histone/Heart

closest-features --closest  --delim '\t' --dist ../../gencode/gencode.v19.annotation.ss ./bed/H3K27ac/H3K27ac_summits.bed |awk '($NF<3000&&$NF>-3000){print $NF}' -  >./enh_find/27ac.dist
closest-features --closest  --delim '\t' --dist ../../gencode/gencode.v19.annotation.ss ./bed/H3K4me1/H3K4me1_summits.bed |awk '($NF<3000&&$NF>-3000){print $NF}' -  >./enh_find/me1.dist
closest-features --closest  --delim '\t' --dist ../../gencode/gencode.v19.annotation.ss ./bed/H3K4me3/H3K4me3_summits.bed |awk '($NF<3000&&$NF>-3000){print $NF}' -  >./enh_find/me3.dist
closest-features --closest  --delim '\t' --dist ../../gencode/gencode.v19.annotation.ss ../../DNase/Heart/Dnase_ENCODE.bed |awk '($NF<3000&&$NF>-3000){print $NF}' -  >./enh_find/dnase.dist


cd /media/ding/000B49000006264C/eRNA_project/histone/Heart
bwtool agg 3000:3000 ../../gencode/protein_lncRNA_sorted.ss ./bg/H3K4me1.bw,./bg/H3K4me3.bw,./bg/H3K27ac.bw,../../DNase/Heart/Dnase_ENCODE.bigwig -firstbase  /dev/stdout > ./enh_find/result.mean_signal





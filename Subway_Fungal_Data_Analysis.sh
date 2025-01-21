[TOC]
#!/bin/bash
# 作者 Xin Zhou
#日期2024.08.24
#使用循环把不同文件夹中文件批量移出来的linux代码
mkdir -p seq
filename=`ls`
for i in  “${filename[@]}”
do
mv "${i}"/*.gz seq/
done

#美格基因数据-2022最新对文件夹中文件批量改名
mkdir rename
for i in `/bin/ls *.gz`; do j=`echo $i | cut -d '_' -f 1,4-5`; cp $i rename/$j; done
cd rename
rename 's/_1//' *.gz

db=../../Public_Scripts

# ${db}
#2. 样本信息，即实验设计 metadata.tsv，保存在最终结果result目录
mkdir -p result
#3. 分析流程pipeline.sh，每个项目复制一份，再进行个性化修改
#4. 创建临时文件存储目录，分析结束可删除
mkdir -p temp
 ### 1.1. metadata.tsv 实验设计文件
    #cat查看前3行，-A显示符号
cat -A result/metadata.txt | head -n3
     #windows用户如果结尾有^M，运行sed命令去除，并cat -A检查结果
 sed -i 's/\r/\n/' result/metadata.txt
cat -A result/metadata.txt | head -n3
### 1.2. seq/*.fq.gz 原始测序数据
    # less按页查看，空格翻页、q退出；head查看前10行，-n指定行
gunzip seq/*.gz


for i in `tail -n+2 result/metadata.txt | cut -f 1`;do
  usearch11 -fastq_mergepairs seq/${i}.R1.fq -reverse seq/${i}.R2.fq \
  -fastqout temp/${i}.merged.fq -relabel ${i}. 
done

# 合并所有样品至同一文件
cat temp/*.merged.fq > temp/all.fq
ls -l temp/all.fq
# 删除中间文件
rm temp/*.merged.fq
# 压缩原始文件节省空间
gzip seq/*

# 3. 切除引物与质控 Cut primers and quality filter
# Cut barcode 10bp + V5 19bp in left and V7 18bp in right
usearch11 -fastx_truncate temp/all.fq \
  -stripleft 10 -stripright 10 \
  -fastqout temp/stripped.fq

# 4. 去冗余与生成OTUs Dereplication and cluster otus
# 质量控制fastq filter, keep reads error rates less than 1%
nohup usearch11 -fastq_filter temp/all.fq \
  -fastq_maxee_rate 0.01 \
  -fastaout temp/filtered.fa &

### 4.1 序列去冗余 Dereplication
    # 并添加miniuniqusize最小为10或1/1M，去除低丰度噪音并增加计算速度
# 4. 去冗余与生成OTUs Dereplication and cluster otus
nohup usearch11 -fastx_uniques  temp/filtered.fa \
-fastaout temp/uniques.fa -minuniquesize 16 -sizeout -relabel Uni & 
### 4.2 聚类OTU/去噪ASV Cluster OTUs / denoise ASV
nohup usearch11 -unoise3 temp/uniques.fa \
  -zotus temp/zotus.fa &
  
 ################################################################################
#修改序列名：Zotu为改为ASV方便识别
sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa
head -n 2 temp/otus.fa
 ################################################################################
#修改序列名：Zotu为改为ASV方便识别
# Win vsearch结果添加了windows换行符^M需删除，mac不要执行此命令
sed -i 's/\r//g' temp/otus.fa

mkdir -p result/raw
# cp temp/otus.fa result/raw/otus.fa
usearch11 -uchime2_ref temp/otus.fa -db ${db}/usearch/UNITE_modified_new2021.fasta \
     -chimeras temp/otus_chimeras.fa -strand plus -mode high_confidence -threads 30
#00:07 170Mb   100.0% Chimeras 45/2920 (1.5%), in db 248 (8.5%), not matched 2627 (90.0%)
#获得非嵌合体序列ID
cat temp/otus.fa temp/otus_chimeras.fa| grep '>' | sort | uniq -u \
  | sed 's/>//' > temp/non_chimeras.id
#筛选非嵌合体，2875 found
usearch11 -fastx_getseqs temp/otus.fa -labels temp/non_chimeras.id \
    -fastaout result/raw/otus.fa
# Win vsearch结果添加了windows换行符^M需删除，mac不要执行此命令
 sed -i 's/\r//g' result/raw/otus.fa 
 ### 5.1 生成特征表 Creat Feature table
nohup usearch11 -otutab temp/filtered.fa -otus result/raw/otus.fa \
   -otutabout result/raw/otutab.txt -threads 20 &
### 5.2 物种注释-去除质体和非细菌/古菌并统计比例(可选) Remove plastid and non-Bacteria

nohup usearch11 -sintax result/raw/otus.fa -db ${db}/usearch/SILVA_138/SILVA_16S_V138.fasta \
-strand both -tabbedout result/raw/otus.sintax -sintax_cutoff 0.2 &

mv result/otutab.txt result/raw
##基于NCBI进行物种注释
#/mnt/data/hanyunpeng
#blastn -max_target_seqs 10 -db /data/db/NCBI_nt/nt -out aspergillus_likeout -query aspergillus_like.fas -num_threads 30 outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids salltitles"
#blast_formatter -archive "othersout" -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids salltitles" > othersout1
#blastn -max_target_seqs 10 -db /data/db/NCBI_nt/nt -out others2out -query others2.fas -num_threads 30 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids salltitles" > othersout
blastn -max_target_seqs 10 -db /data/db/NCBI_nt/nt -out penicillium5out -query penicillium5.fas -num_threads 30 -outfmt "7 qseqid pident salltitles"

nohup blastn -max_target_seqs 5 -db /data/db/NCBI_nt/nt -out Wuwei_Fungi_NCBI.blast \
-query otus.fa -num_threads 20 \
-outfmt "6 qseqid qlen qstart qend sseqid stitle slen sstart send qcovs bitscore score evalue pident" &

grep -v '^#' penicillium5out | awk '{print$1,$2,$3,$4}' > penicilliumout5 && rm -f penicillium5out
# 原始特征表行数,f:,g:
wc -l result/otutab.txt
#R脚本选择细菌古菌(真核)、去除叶绿体、线粒体并统计比例；输出筛选并排序的OTU表
#输入为OTU表result/raw/otutab.txt和物种注释result/raw/otus.sintax
#输出筛选并排序的特征表result/otutab.txt和
#统计污染比例文件result/raw/otutab_nonBac.txt和过滤细节otus.sintax.discard
#将来更新为按taxonomy自由过滤的参数，适合细菌、真菌
#Rscript ${db}/script/otutab_filter_nonFungi.R
 #将来更新为按taxonomy自由过滤的参数，适合细菌、真菌
db=../../Public_Scripts
##过滤真菌筛选完还有 (5,860,384 reads)
Rscript ${db}/script/otutab_filter_nonFungi.R -h
Rscript ${db}/script/otutab_filter_nonFungi.R \
      --input result/raw/otutab.txt \
      --taxonomy result/raw/otus.sintax \
      --output result/otutab.txt\
      --stat result/raw/otutab_nonFunc.stat \
      --discard result/raw/otus.sintax.discard
      
# 筛选后特征表行数
wc -l result/otutab.txt
usearch11 -otutab_trim result/otutab.txt -min_otu_size 16 -output result/otutab1.txt
mv result/otutab1.txt result/otutab.txt
#提取ID用于提取序列
cut -f 1 result/otutab.txt | sed '1 s/#OTU ID/OTUID/' > result/otu.id
wc -l result/otu.id
#筛选对应OTU序列
usearch11 -fastx_getseqs result/raw/otus.fa -labels result/otu.id \
-fastaout result/otus.fa

#过滤特征表对应序列注释
awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}'\
    result/raw/otus.sintax result/otutab.id \
    > result/otus.sintax
#补齐末尾列
sed -i 's/\t$/\td:Unassigned/' result/otus.sintax

## 8. 物种注释结果分类汇总 Taxonomy summary
#OTU对应物种注释2列格式：去除sintax中置信值，只保留物种注释，替换:为_，删除引号
cut -f 1,4 result/otus.sintax \
  |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' \
   > result/taxonomy2.txt
head -n13 result/taxonomy2.txt

#OTU对应物种8列格式：注意注释是非整齐
#生成物种表格OTU/ASV中空白补齐为Unassigned
awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
  split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
  print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
result/taxonomy2.txt > temp/otus.tax
awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
  split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
  print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
result/taxonomy2.txt > temp/otus.tax
sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
> result/taxonomy.txt
head -n3 result/taxonomy.txt

grep ">" result/otus.fa | wc -l
#可选统计方法：OTU表简单统计 Summary OTUs table
usearch11 -otutab_stats result/otutab.txt \
  -output result/otutab.stat
cat result/otutab.stat
#注意最小值、分位数，或查看result/raw/otutab_nonBac.stat中样本详细数据量，用于重采样
### 5.3 等量抽样标准化 normlize by subsample
#使用vegan包进行等量重抽样，输入reads count格式Feature表result/otutab.txt
#49548170  Reads (49.5M)  466  Samples 31496  OTUs 细菌的抽平数为90531
db=../../Public_Scripts
mkdir -p result/alpha
Rscript ${db}/script/otutab_rare.R --input result/otutab.txt \
  --depth 40907 --seed 1 \
  --normalize result/otutab_rare.txt \
  --output result/alpha/vegan.txt

usearch11 -otutab_stats result/otutab_rare.txt \
      -output result/otutab_rare.stat
cat result/otutab_rare.stat  
#统计门纲目科属，使用 rank参数 p c o f g，为phylum, class, order, family, genus缩写
mkdir -p result/tax
for i in p c o f g ;do
  usearch11 -sintax_summary result/otus.sintax \
  -otutabin result/otutab_rare.txt -rank ${i} \
  -output result/tax/sum_${i}.txt
done
sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' result/tax/sum_*.txt
# 列出所有文件
ls -sh result/tax/sum_*.txt
head -n3 result/tax/sum_g.txt

##把物种注释信息分号替换成TAB分隔
cat result/taxonomy2.txt | tr -s ";" "\t" > result/rep_seqs_tax.txt
head result/rep_seqs_tax.txt

#计算各特征的均值，有组再求分组均值，需根据实验设计metadata.txt修改组列名
#删除文件最后一列
#awk '{print $NF}' result/tax/test.txt
awk 'BEGIN{IFS='\t'}{$NF="";print $0}' result/tax/sum_g.txt | head -n 12 > result/tax/abundant_genus.txt
##删除每行最后的空格
sed -i 's/[ ]*$//g' result/tax/abundant_genus.txt
##删除第二行未注释文件
#sed -i '2d' result/tax/abundant_genus.txt
cat result/tax/abundant_genus.txt | tr -s ' ' '\t' > result/tax/abundant_genus1.txt
cat result/tax/abundant_genus1.txt
mv result/tax/abundant_genus1.txt result/tax/abundant_genus.txt
cat result/tax/abundant_genus.txt
##删除文件中All 这列
cat result/tax/final_g_group.txt | awk -F '\t' '{$2=null;print $0}' | tr -s ' ' '\t' > result/tax/final_g_group1.txt
mv result/tax/final_g_group1.txt result/tax/final_g_group.txt

##提取指定分组中的OTU表
mkdir -p Intercity/result
grep "Intercity" result/metadata.txt | cut -f 1 > Intercity/result/SampleID.txt

usearch11 -otutab_sample_subset result/otutab.txt \
-labels Intercity/result/SampleID.txt -output Intercity/result/otutab1.txt

#usearch参数解读：-min_sample_size reads于5000的；min_count；-min_freq 和所有样品比；-min_otu_size 每个OTU
#-min_otu_freq Minimum size for an OTU as a fraction of all OTUs
usearch11 -otutab_trim Intercity/result/otutab1.txt -min_otu_size 2 -output Intercity/result/otutab.txt
#可选统计方法：OTU表简单统计 Summary OTUs table
usearch11 -otutab_stats Intercity/result/otutab.txt \
   -output Intercity/result/otutab.stat
cat Intercity/result/otutab.stat
#提取ID用于提取序列
cut -f 1 Intercity/result/otutab.txt | sed '1 s/#OTU ID/OTUID/' > Intercity/result/otu.id
wc -l Intercity/result/otu.id
#筛选对应OTU序列
usearch11 -fastx_getseqs result/raw/otus.fa -labels Intercity/result/otu.id \
-fastaout Intercity/result/otus.fa
mkdir -p Intercity/result/raw/
## 筛选OTU对物种注释
awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/otus.sintax \
Intercity/result/otu.id > Intercity/result/raw/otus.sintax
head Intercity/result/raw/otus.sintax
awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/taxonomy.txt \
Intercity/result/otu.id > Intercity/result/taxonomy.txt
cp Intercity/result/raw/otus.sintax Intercity/result/otus.sintax

### 6.1. 计算多样性指数 Calculate alpha diversity index
usearch11 -alpha_div result/otutab_rare.txt \
  -output result/alpha/alpha.txt
### 6.2. 计算稀释过程的丰富度变化 Rarefaction
#Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)-method fast  https://drive5.com/usearch/manual/cmd_otutab_subsample.html
usearch10 -alpha_div_rare result/otutab.txt \
  -output result/alpha/alpha_rare.txt -method without_replacement
##Beta多样性 Beta diversity
mkdir -p result/beta/
#生成5种距离矩阵：bray_curtis, euclidean, jaccard, manhatten, unifrac, 3s
usearch11 -beta_div result/otutab_rare.txt -filename_prefix result/beta/

##核心微生物分析
usearch11 -calc_distmx result/otus.fa -tabbedout result/distmx.txt -maxdist 0.8 -termdist 0.9
  
usearch11 -otutab_core result/otutab.txt -distmxin result/distmx.txt \
  -sintaxin result/otus.sintax -tabbedout result/core.txt
### 6.3. 筛选各组高丰度菌用于比较
#输入文件为feautre表result/otutab.txt，实验设计metadata.txt
#输出为特征表按组的均值-一个实验可能有多种分组方式
Rscript ${db}/script/otu_mean.R --input result/otutab_rare.txt \
   --design result/metadata.txt \
   --group Seasonal --thre 0 \
   --output result/otutab_mean_Seasonal.txt
head -n3 result/otutab_mean_Seasonal.txt
#######################################################################
#如以平均丰度频率高于千分之一(0.1%)为筛选标准，得到每个组的OTU组合
awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
    else {for(i=2;i<=NF;i++) if($i>0.0001) print $1, a[i];}}' \
    result/otutab_mean_Seasonal.txt > result/alpha/otutab_mean_Seasonal.txt
head result/alpha/otutab_mean_Seasonal.txt

db=../../../../Public_Scripts
### 1.2 稀释曲线
Rscript ${db}/script/alpha_rare_curve.R \
  --input result/alpha/alpha_rare.txt --design result/metadata.txt \
  --group Station --output result/alpha/ \
  --width 180 --height 120
# 23、R语言多样性和物种分析
## 1. Alpha多样性
Rscript ${db}/script/alpha_boxplot.R -h
# 完整参数，多样性指数可选richness chao1 ACE shannon simpson invsimpson
Rscript ${db}/script/alpha_boxplot.R --alpha_index richness \
  --input result/alpha/vegan.txt --design result/metadata.txt \
  --group POSITION --output result/alpha/ \
  --width 90 --height 60
# 使用循环绘制6种常用指数
for i in `head -n1 result/alpha/vegan.txt|cut -f 2-`;do
  Rscript ${db}/script/alpha_boxplot.R --alpha_index ${i} \
     --input result/alpha/vegan.txt --design result/metadata.txt \
     --group POSITION --output result/alpha/ \
    --width 180 --height 120
done
### 1.3 多样性维恩图
# 四组比较，图和代码见输入文件目录，运行目录为当前项目根目录
bash ${db}/script/sp_vennDiagram.sh \
  -f result/alpha/otu_group_exist_Station_Type.txt \
  -a Intercity  -b Intercity -c Intercity \
  -w 3 -u 3 \
  -p Station_Venn

## 2. Beta多样性
### 2.1 距离矩阵热图pheatmap
# 以bray_curtis为例，-f输入文件,-h是否聚类TRUE/FALSE,-u/v为宽高英寸
bash ${db}/script/sp_pheatmap.sh \
  -f result/beta/bray_curtis.txt \
  -H 'TRUE' -u 16 -v 12
# 添加分组注释，如2，4列的基因型和地点
cut -f 1,3 result/metadata.txt > temp/group.txt
# -P添加行注释文件，-Q添加列注释
bash ${db}/script/sp_pheatmap.sh \
  -f result/beta/bray_curtis.txt \
  -H 'TRUE' -u 16 -v 12 \
  -P temp/group.txt -Q temp/group.txt
# 距离矩阵与相关类似，可尝试corrplot或ggcorrplot绘制更多样式
# - [绘图相关系数矩阵corrplot](http://mp.weixin.qq.com/s/H4_2_vb2w_njxPziDzV4HQ)
# - [相关矩阵可视化ggcorrplot](http://mp.weixin.qq.com/s/AEfPqWO3S0mRnDZ_Ws9fnw)

## 3. 物种组成Taxonomy
### 3.1 堆叠柱状图Stackplot
# 以门(p)水平为例，结果包括output.sample/group.pdf两个文件
Rscript ${db}/script/tax_stackplot.R \
  --input result/tax/sum_${i}.txt --design result/metadata.txt \
  --group group --output result/tax/sum_${i}.stackplot \
  --legend 2 --width 180 --height 120
# 批量绘制输入包括p/c/o/f/g共5级
for i in p c o f g s; do
Rscript ${db}/script/tax_stackplot.R \
  --input result/tax/sum_${i}.txt --design result/metadata.txt \
  --group POSITION --output result/tax/sum_${i}.stackplot \
  --legend 13 --width 180 --height 120; 
done

### 3.2 弦图(圈图)circlize
# 以纲(class,c)为例，绘制前5组
i=o
Rscript ${db}/script/tax_circlize.R \
  --input result/tax/sum_${i}.txt --design result/metadata.txt \
  --group group --legend 10
# 结果位于当前目录circlize.pdf(随机颜色)，circlize_legend.pdf(指定颜色+图例)
# 移动并改名与分类级一致
mv circlize.pdf result/tax/sum_${i}.circlize.pdf
mv circlize_legend.pdf result/tax/sum_${i}.circlize_legend.pdf

# 24、差异比较
## 24-1. R语言差异分析
### 1.1 差异比较Difference comparison
compare="Warm-Cold"
Rscript ${db}/script/compare.R \
  --input result/otutab_rare.txt --design result/metadata.txt \
  --group Seasonal1 --compare ${compare} --threshold 0.001 \
  --method edgeR --pvalue 0.05 --fdr 0.2 \
  --output result/compare/
### 1.2 火山图

# 输入compare.R的结果，输出火山图带数据标签，可指定图片大小
Rscript ${db}/script/compare_volcano.R \
   --input result/compare/${compare}.txt \
   --output result/compare/${compare}.volcano.pdf \
   --width 120 --height 120

### 1.3 热图
# 输入compare.R的结果，筛选列数，指定元数据和分组、物种注释，图大小英寸和字号
bash ${db}/script/compare_heatmap.sh -i result/compare/${compare}.txt -l 7 \
       -d result/metadata.txt -A Seasonal1 \
       -t result/taxonomy.txt \
       -w 10 -h 14 -s 8 \
       -o result/compare/${compare}.txt

### 1.4 曼哈顿图
# i差异比较结果,t物种注释,p图例,w宽,v高,s字号,l图例最大值
bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
   -t result/taxonomy.txt \
   -p result/tax/sum_p.txt \
   -w 180 -v 60 -s 7 -l 10 \
   -o result/compare/${compare}
# 上图只有6个门，切换为纲c和-L Class展示细节
bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
       -t result/taxonomy.txt \
       -p result/tax/sum_c.txt \
       -w 183 -v 80 -s 7 -l 10 -L Class \
       -o result/compare/${compare}.txt

# 显示完整图例，再用AI拼图
bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
       -t result/taxonomy.txt \
       -p result/tax/sum_c.txt \
       -w 183 -v 149 -s 7 -l 10 -L Class \
       -o result/compare/${compare}.manhattan.c.legend.pdf
### 1.5 单个特征的绘制
  
# 差异OTU细节展示
  Rscript ${db}/script/alpha_boxplot.R --alpha_index Lactobacillus \
    --input result/tax/sum_g.txt --design result/metadata.txt \
    --transpose TRUE --scale TRUE \
    --width 89 --height 59 \
    --group group --output result/tax/

# 差属细节展示
    Rscript ${db}/script/alpha_boxplot.R --alpha_index Filobasidium \
      --input result/tax/sum_g.txt --design result/metadata.txt \
      --transpose TRUE \
      --width 89 --height 59 \
      --group group1 --output result/compare/
## 24-2. STAMP输入文件准备
Rscript ${db}/script/format2stamp.R -h
mkdir -p result/stamp
Rscript ${db}/script/format2stamp.R --input result/otutab_rare.txt \
   --taxonomy result/taxonomy.txt --threshold 0.001 \
   --output result/tax

## 24-3. LEfSe输入文件准备
### 3.1. 命令行生成文件
# 可选命令行生成输入文件
Rscript ${db}/script/format2lefse.R -h
Rscript ${db}/script/format2lefse.R --input result/otutab.txt \
  --taxonomy result/taxonomy.txt --design result/metadata.txt \
  --group StationType --threshold 0.0001 \
  --output result/LEfSe_Station_Type

##制作三元图输入文件
#合并物种注释和丰度为注释文件
mkdir -p result/compare/Ternary_plot
cut -f 2- otutab_high.mean > temp

下面就是如何使用sed往一个文件顶部添加一行的方法：

##在taxonomy第一行插入名称
sed '1i OTUId \t Silva_taxonomy' result/taxonomy2.txt >result/taxonomy3.txt
##合并物种注释表和OTU 表
paste result/taxonomy3.txt result/otutab_rare.txt > result/otu_table.txt
head -n 3 result/otu_table.txt
##删除文件中#OTUID这列
cat result/otu_table.txt | awk -F '\t' '{$3=null;print $0}' | tr -s ' ' '\t' > result/otu_table1.txt
##查找替换分类注释符
sed 's/k__/D_0__/g' result/otu_table1.txt > result/otu_table.txt
sed 's/p__/D_1__/g' result/otu_table.txt > result/otu_table1.txt
sed 's/c__/D_2__/g' result/otu_table1.txt > result/otu_table.txt
sed 's/o__/D_3__/g' result/otu_table.txt > result/otu_table1.txt
sed 's/f__/D_4__/g' result/otu_table1.txt > result/otu_table.txt
sed 's/g__/D_5__/g' result/otu_table.txt > result/otu_table1.txt
sed 's/s__/D_6__/g' result/otu_table1.txt > result/otu_table.txt

Rscript ${db}/script/tax_circlize.R \
  --input KEGG.Pathway.raw.txt --design ../metadata.txt \
  --group group1 --legend 11
# 移动并改名与分类级一致
mv circlize.pdf result/sum_KEGG.Pathway.circlize.pdf
mv circlize_legend.pdf result/sum_KEGG.Pathway.circlize_legend.pdf

mv result/otu_table.txt result/Ternary_plot/otu_table.txt


# 33、MachineLearning机器学习
# RandomForest包使用的R代码见33MachineLearning目录中的RF_classification和RF_regression
## Silme2随机森林/Adaboost使用代码见33MachineLearning目录中的slime2，或附录5
#USEARCH forest_train 要求的输入文件这里定义为  feature_table,  
#对传统的 “OTU_table” 进行转置（行为样本，列为OTU）
#metadata.txt为样本分类信息，第一列：样本，第二列：分组
tsv-utils  transpose OTU_table.txt            \
| tsv-utils annotation                        \
       <(tsv-utils  add_headline              \
           "#OTU\tCatalog" metadata.txt)  -   \
 | sed 's/#OTU ID/Samples/'                   \
 > feature_table.txt
 #执行forest_train; forest_train 其实就是对随机森林分类训练参数
 #对参数文件的描述见： https://www.drive5.com/usearch/manual/forest_params.html
 usearch11  -forest_train     \
   feature_table.txt   -trees 100 -forestout forest.txt
 #获取OTU排序信息
#通过提取训练参数文件中的变量（OTU）重要性值，并对其进行排序。
grep -w "^varimp" forest.txt                  \
|cut -f3,4                                    \
|sort -rgk2                                   \
| tsv-utils annotation sintax.txt -           \
| tsv-utils add_headline                      \
     "#OTU\tImportant\tLineage" -             \
>rank.txt 

#输出为特征表按组的均值-一个实验可能有多种分组方式
Rscript ${db}/script/otu_mean.R --input result/Pathogen_otutab_rare.txt \
   --design result/metadata.txt \
   --group Station --thre 0 \
   --output result/otutab_mean_Station.txt
head -n3 result/otutab_mean_Station.txt
#######################################################################
#如以平均丰度频率高于千分之一(0.1%)为筛选标准，得到每个组的OTU组合
awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
    else {for(i=2;i<=NF;i++) if($i>0.0001) print $1, a[i];}}' \
    result/otutab_mean_Seasonal.txt > result/alpha/otutab_mean_Seasonal.txt
head result/alpha/otutab_mean_Seasonal.txt

# 起始文件为 result/tree目录中 otus.fa(序列)、annotation.txt(物种和相对丰度)文件
cd /mnt/sdd/ZhouXin/Tomato_project/DataProcess4Publication/Field_Bacteria_16S/Bulk
mkdir -p result/OTU_tree
# Muscle软件进行序列对齐
nohup muscle -in result/otus.fa -out result/OTU_tree/otus_aligned.fas &
##必须增加一步，修剪序列
#trimAL 软件进行低质量以及高变异度的序列的过滤和修剪，代码如下：
trimal -in result/OTU_tree/otus_aligned.fas -out result/OTU_tree/otus_aligned_trimed.fa -gt 0.95
### 利用IQ-TREE快速构建ML进化树
nohup iqtree -s result/OTU_tree/otus_aligned_trimed.fa \
        -bb 5000 -redo -alrt 5000 -nt AUTO -T 40  \
        -pre result/OTU_tree/otus &
        
####Funguild分析

cut -f 1 result/growthForm_sintax.txt | sed '1 s/#OTU ID/OTUID/' > result/growthForm.id

#根据ASVs ID筛选子OTU表
usearch11 -otutab_otu_subset result/otutab.txt -labels result/growthForm.id \
-output result/otutab_growthForm.txt

usearch11 -sintax_summary result/growthForm_sintax.txt \
  -otutabin result/otutab_growthForm.txt -rank g \
  -output result/sum_growthForm.txt

cut -f 1 result/guild_sintax.txt | sed '1 s/#OTU ID/OTUID/' > result/guild.id
#根据guild_sintax ID筛选子OTU表
usearch11 -otutab_otu_subset result/otutab.txt -labels result/guild.id \
-output result/otutab_guild.txt

usearch11 -sintax_summary result/guild_sintax.txt \
  -otutabin result/otutab_guild.txt -rank g \
  -output result/sum_guild.txt


cut -f 1 result/trait_sintax.txt | sed '1 s/#OTU ID/OTUID/' > result/trait.id
#根据guild_sintax ID筛选子OTU表
usearch11 -otutab_otu_subset result/otutab.txt -labels result/trait.id \
-output result/otutab_trait.txt

usearch11 -sintax_summary result/trait_sintax.txt \
  -otutabin result/otutab_trait.txt -rank g \
  -output result/sum_trait.txt

db=../../../../Public_Scripts
##LEfSe输入文件准备
Rscript ${db}/script/format2lefse.R -h
Rscript ${db}/script/format2lefse.R --input result/otutab_guild.txt \
  --taxonomy Lefse/guild_taxonomy.txt --design result/metadata.txt \
  --group Station_Type --threshold 0.0001 \
  --output Lefse/Station_Type


# 差异OTU细节展示
  Rscript ${db}/script/alpha_boxplot.R --alpha_index SoftRot \
    --input result/sum_trait.txt --design result/metadata.txt \
    --transpose TRUE --scale TRUE \
    --width 89 --height 59 \
    --group group --output result/  


#############################################################
#过滤特征表对应序列注释
awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}'\
    ../otus.sintax ./Core_OTUID.txt \
    > ./otus.sintax
#补齐末尾列
sed -i 's/\t$/\td:Unassigned/' ./otus.sintax

## 8. 物种注释结果分类汇总 Taxonomy summary
#OTU对应物种注释2列格式：去除sintax中置信值，只保留物种注释，替换:为_，删除引号
cut -f 1,4 ./otus.sintax \
  |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' \
   > ./taxonomy2.txt
head -n13 ./taxonomy2.txt

#OTU对应物种8列格式：注意注释是非整齐
#生成物种表格OTU/ASV中空白补齐为Unassigned
awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
  split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
  print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
./taxonomy2.txt > ./otus.tax
awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
  split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
  print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
./taxonomy2.txt > ./otus.tax
sed 's/;/\t/g;s/.__//g;' ./otus.tax|cut -f 1-8 | \
sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
> ./taxonomy.txt
head -n3 ./taxonomy.txt

#统计门纲目科属，使用 rank参数 p c o f g，为phylum, class, order, family, genus缩写
mkdir -p ./tax
for i in p c o f g ;do
  usearch11 -sintax_summary ./otus.sintax \
  -otutabin ./core_otutab.txt -rank ${i} \
  -output ./tax/sum_${i}.txt
done
sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' ./tax/sum_*.txt
# 列出所有文件
ls -sh ./tax/sum_*.txt
head -n3 ./tax/sum_g.txt

##把物种注释信息分号替换成TAB分隔
cat ./taxonomy2.txt | tr -s ";" "\t" > ./rep_seqs_tax.txt
head ./rep_seqs_tax.txt

#计算各特征的均值，有组再求分组均值，需根据实验设计metadata.txt修改组列名
#删除文件最后一列
#awk '{print $NF}' result/tax/test.txt
awk 'BEGIN{IFS='\t'}{$NF="";print $0}'  ./tax/sum_g.txt | head -n 15 > ./tax/abundant_genus.txt
##删除每行最后的空格
sed -i 's/[ ]*$//g' ./tax/abundant_genus.txt
##删除第二行未注释文件
#sed -i '2d' result/tax/abundant_genus.txt
cat ./tax/abundant_genus.txt | tr -s ' ' '\t' > ./tax/abundant_genus1.txt
cat ./tax/abundant_genus1.txt
mv ./tax/abundant_genus1.txt ./tax/abundant_genus.txt
cat ./tax/abundant_genus.txt
###############################################################################
#################################核心物种分析
## 2. 构建进化树

mkdir -p result/Core_OTU_Tree
cd result/Core_OTU_Tree

#根据ASVs ID筛选子OTU表
usearch11 -otutab_otu_subset ../raw/otutab.txt -labels ./Core_OTUID.txt -output ./result/core_otutab.txt
#筛选指定差异菌对应OTU序列
usearch11 -fastx_getseqs ../otus.fa -labels ./Core_OTUID.txt \
    -fastaout ./Core_otus.fa
head -n 2 ./Core_otus.fa
grep ">" ./Core_otus.fa |wc -l
## 筛选OTU对物种注释
awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../taxonomy.txt \
./Core_OTUID.txt > ./Core_otutab_high.tax

#获得OTU对应组均值，用于样本热图
#依赖之前otu_mean.R计算过按Group分组的均值
awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../otutab_mean_Station_Type.txt ./Core_OTUID.txt \
    | sed 's/#OTU ID/OTUID/' > ./otutab_high_Station.mean


awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../otutab_mean_Seasonal.txt ./Core_OTUID.txt \
    | sed 's/#OTU ID/OTUID/' > ./otutab_mean_Seasonal.mean
#sed -i '2d' ./otutab_high.mean
head -n3 ./otutab_high_Station.mean
#合并物种注释和丰度为注释文件
cut -f 2- ./otutab_high_Station.mean > temp1
cut -f 2- ./otutab_mean_Seasonal.mean > temp2
paste ./Core_otutab_high.tax temp1 temp2> ./annotation1.txt
head -n 3 ./annotation1.txt
##添加另外一个分组信息
cut -f 2- ./otutab_high_Station.mean > temp1
paste ./annotation1.txt temp1 > ./annotation2.txt


# 起始文件为 result/tree目录中 otus.fa(序列)、annotation.txt(物种和相对丰度)文件
# Muscle软件进行序列对齐，3s
time muscle -in Core_otus.fa -out otus_aligned.fas
#trimAL 软件进行低质量以及高变异度的序列的过滤和修剪，代码如下：
trimal -in otus_aligned.fas -out otus_aligned_trimed.fa -gt 0.95
### 方法1. 利用IQ-TREE快速构建ML进化树，2m
mkdir -p iqtree
iqtree2 -s  otus_aligned_trimed.fa \
        -bb 1000 -redo -alrt 1000 -nt AUTO \
        -pre iqtree/Core_otus
## 3. 进化树美化
# 访问http://itol.embl.de/，上传otus.nwk，再拖拽下方生成的注释方案于树上即美化
## 方案1. 外圈颜色、形状分类和丰度方案
# annotation.txt OTU对应物种注释和丰度，
# -a 找不到输入列将终止运行（默认不执行）-c 将整数列转换为factor或具有小数点的数字，-t 偏离提示标签时转换ID列，-w 颜色带，区域宽度等， -D输出目录，-i OTU列名，-l OTU显示名称如种/属/科名，
db=../../../../Public_Scripts
Rscript ${db}/script/table2itol.R -a -c double -D plan1 -i OTUID -l Family -t %s -w 0.5 annotation.txt
# 生成注释文件中每列为单独一个文件
## 方案2. 生成丰度柱形图注释文件
Rscript ${db}/script/table2itol.R -a -d -c none -D plan2 -b Phylum -i OTUID -l Family -t %s -w 0.5 annotation.txt
## 方案3. 生成热图注释文件
Rscript ${db}/script/table2itol.R -c keep -D plan3 -i OTUID -t %s otutab_mean_Seasonal.mean
## 方案4. 将整数转化成因子生成注释文件
Rscript ${db}/script/table2itol.R -a -c factor -D plan4 -i OTUID -l Genus -t %s -w 0 annotation.txt
# 返回工作目录
#可选统计方法：OTU表简单统计 Summary OTUs table
usearch11 -otutab_stats result/core_otutab.txt \
  -output result/otutab.stat
cat result/otutab.stat

db=../../../../Public_Scripts
mkdir -p result/alpha
Rscript ${db}/script/otutab_rare.R --input result/core_otutab.txt \
  --depth 25060 --seed 1 \
  --normalize result/otutab_rare.txt \
  --output result/alpha/vegan.txt


db=../../../../Public_Scripts
##LEfSe输入文件准备
Rscript ${db}/script/format2lefse.R -h
Rscript ${db}/script/format2lefse.R --input result/otutab_rare.txt \
  --taxonomy result/taxonomy.txt --design ../metadata.txt \
  --group Seasonal --threshold 0.0001 \
  --output Lefse/Seasonal

##制作三元图输入文件
#合并物种注释和丰度为注释文件
mkdir -p result/compare/Ternary_plot
cut -f 2- otutab_high.mean > temp
下面就是如何使用sed往一个文件顶部添加一行的方法：
##在taxonomy第一行插入名称
sed '1i OTUId \t Silva_taxonomy' result/taxonomy2.txt >result/taxonomy3.txt
##合并物种注释表和OTU 表
paste result/taxonomy3.txt result/otutab_rare.txt > result/otu_table.txt
head -n 3 result/otu_table.txt
##删除文件中#OTUID这列
cat result/otu_table.txt | awk -F '\t' '{$3=null;print $0}' | tr -s ' ' '\t' > result/otu_table1.txt
##查找替换分类注释符
sed 's/k__/D_0__/g' result/otu_table1.txt > result/otu_table.txt
sed 's/p__/D_1__/g' result/otu_table.txt > result/otu_table1.txt
sed 's/c__/D_2__/g' result/otu_table1.txt > result/otu_table.txt
sed 's/o__/D_3__/g' result/otu_table.txt > result/otu_table1.txt
sed 's/f__/D_4__/g' result/otu_table1.txt > result/otu_table.txt
sed 's/g__/D_5__/g' result/otu_table.txt > result/otu_table1.txt
sed 's/s__/D_6__/g' result/otu_table1.txt > result/otu_table.txt

##提取指定分组中的OTU表，Tomato
mkdir -p Intercity/result
grep "Intercity" ../metadata.txt | cut -f 1 > Intercity/result/SampleID.txt

usearch11 -otutab_sample_subset result/core_otutab.txt \
-labels Intercity/result/SampleID.txt -output Intercity/result/otutab1.txt

#usearch参数解读：-min_sample_size reads于5000的；min_count；-min_freq 和所有样品比；-min_otu_size 每个OTU
#-min_otu_freq Minimum size for an OTU as a fraction of all OTUs
usearch11 -otutab_trim Intercity/result/otutab1.txt -min_otu_size 2 -output Intercity/result/otutab.txt

#提取ID用于提取序列
cut -f 1 Intercity/result/otutab.txt | sed '1 s/#OTU ID/OTUID/' > Intercity/result/otu.id
wc -l Intercity/result/otu.id
#筛选对应OTU序列
usearch11 -fastx_getseqs ../raw/otus.fa -labels Intercity/result/otu.id \
-fastaout Intercity/result/otus.fa
mkdir -p Intercity/result/raw/
## 筛选OTU对物种注释
awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../otus.sintax \
Intercity/result/otu.id > Intercity/result/raw/otus.sintax
head Intercity/result/raw/otus.sintax
awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../taxonomy.txt \
Intercity/result/otu.id > Intercity/result/taxonomy.txt
cp Intercity/result/raw/otus.sintax Intercity/result/otus.sintax

#可选统计方法：OTU表简单统计 Summary OTUs table
usearch11 -otutab_stats Intercity/result/otutab.txt \
   -output Intercity/result/otutab.stat
cat Intercity/result/otutab.stat
#统计门纲目科属，使用 rank参数 p c o f g，为phylum, class, order, family, genus缩写
db=../../../../Public_Scripts
mkdir -p Intercity/result/alpha
Rscript ${db}/script/otutab_rare.R --input Intercity/result/otutab.txt \
  --depth 23745 --seed 1 \
  --normalize Intercity/result/otutab_rare.txt \
  --output Intercity/result/alpha/vegan.txt
  

mkdir -p Intercity/result/tax
for i in p c o f g ;do
  usearch11 -sintax_summary Intercity/result/otus.sintax \
  -otutabin Intercity/result/otutab_rare.txt -rank ${i} \
  -output Intercity/result/tax/sum_${i}.txt
done
sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' Intercity/result/tax/sum_*.txt

###病原真菌分析
db=../../../../Public_Scripts
mkdir -p result/alpha
Rscript ${db}/script/otutab_rare.R --input ./Pathogen_otutab_rare1.txt \
  --depth 5061 --seed 1 \
  --normalize result/Pathogen_otutab_rare.txt \
  --output result/alpha/vegan.txt

usearch11 -otutab_stats result/Pathogen_otutab_rare.txt \
      -output result/otutab_rare.stat
cat result/otutab_rare.stat  

db=../../../../Public_Scripts
##LEfSe输入文件准备
mkdir -p Lefse
Rscript ${db}/script/format2lefse.R -h
Rscript ${db}/script/format2lefse.R --input result/Pathogen_otutab_rare.txt \
  --taxonomy ./taxonomy.txt --design result/metadata.txt \
  --group Station_Type --threshold 0.0001 \
  --output Lefse/Station_Type

grep "Spring" result/metadata.txt > result/metadata_Spring.txt
##添加另外一个分组信息
sed -n '1p' result/metadata.txt > temp1
cat temp1 result/metadata_Spring.txt > temp_combined && mv temp_combined result/metadata_Spring1.txt
mv result/metadata_Spring1.txt result/metadata_Spring.txt

grep "Summer" result/metadata.txt > result/metadata_Summer.txt
##添加另外一个分组信息
sed -n '1p' result/metadata.txt > temp1
cat temp1 result/metadata_Summer.txt > temp_combined && mv temp_combined result/metadata_Summer.txt

grep "Autumn" result/metadata.txt > result/metadata_Autumn.txt
##添加另外一个分组信息
sed -n '1p' result/metadata.txt > temp1
cat temp1 result/metadata_Autumn.txt > temp_combined && mv temp_combined result/metadata_Autumn.txt

grep "Winter" result/metadata.txt > result/metadata_Winter.txt
##添加另外一个分组信息
sed -n '1p' result/metadata.txt > temp1

cat temp1 result/metadata_Winter.txt > temp_combined && mv temp_combined result/metadata_Winter.txt


#输出为特征表按组的均值-一个实验可能有多种分组方式
Rscript ${db}/script/otu_mean.R --input result/Pathogen_otutab_rare.txt \
   --design result/metadata.txt \
   --group Station --thre 0 \
   --output result/otutab_mean_Station.txt
head -n3 result/otutab_mean_Station.txt


## 筛选OTU对物种注释
awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/metadata.txt \
result/AutumnID.txt > result/metadata_Autumn.txt
cat temp1 result/metadata_Autumn.txt > temp_combined && mv temp_combined result/metadata_Autumn1.txt
mv result/metadata_Autumn1.txt result/metadata_Autumn.txt

grep "Winter" result/metadata.txt | cut -f 1 > result/SampleID.txt
usearch11 -otutab_sample_subset result/Pathogen_otutab_rare.txt \
-labels result/SampleID.txt -output result/Pathogen_otutab_Winter.txt
#输出为特征表按组的均值-一个实验可能有多种分组方式
Rscript ${db}/script/otu_mean.R --input result/Pathogen_otutab_Winter.txt \
   --design result/metadata_Winter.txt \
   --group Station --thre 0 \
   --output result/otutab_mean_Station_Winter.txt
head -n3 result/otutab_mean_Station_Winter.txt
#如以平均丰度频率高于千分之一(0.1%)为筛选标准，得到每个组的OTU组合
awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
    else {for(i=2;i<=NF;i++) if($i>0.0001) print $1, a[i];}}' \
    result/otutab_mean_Station.txt > result/alpha/otu_group_exist_Station.txt
head result/alpha/otu_group_exist_Station.txt


awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
    else {for(i=2;i<=NF;i++) if($i>0.0001) print $1, a[i];}}' \
    result/otutab_mean_Station_Winter.txt > result/alpha/otu_group_exist_Station_Winter.txt
head result/alpha/otu_group_exist_Station_Winter.txt
#如以平均丰度频率高于千分之一(0.1%)为筛选标准，得到每个组的OTU组合
awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
    else {for(i=2;i<=NF;i++) if($i>0.0001) print $1, a[i];}}' \
    result/otutab_mean_Seasonal.txt > result/alpha/otu_group_exist_Seasonal.txt
head result/alpha/otu_group_exist_Seasonal.txt
### 1.3 多样性维恩图
# 四组比较，图和代码见输入文件目录，运行目录为当前项目根目录
bash ${db}/script/sp_vennDiagram.sh \
  -f result/alpha/otu_group_exist_Seasonal.txt \
  -a Spring -b Summer -c Autumn -d Winter \
  -w 5 -u 5 \
  -p Venn_Seasonal

bash ${db}/script/sp_vennDiagram.sh \
  -f result/alpha/otu_group_exist_Station_Winter.txt \
  -a Residential -b Intercity -c UrbanHub \
  -w 5 -u 5 \
  -p Venn_Station_Winter   

  
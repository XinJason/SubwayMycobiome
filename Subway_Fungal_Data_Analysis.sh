[TOC]
#!/bin/bash
# 作者 Xin Zhou

mkdir rename
for i in `/bin/ls *.gz`; do j=`echo $i | cut -d '_' -f 1,4-5`; cp $i rename/$j; done
cd rename
rename 's/_1//' *.gz

db=../../Public_Scripts

# ${db}

mkdir -p result

mkdir -p temp


 sed -i 's/\r/\n/' result/metadata.txt
cat -A result/metadata.txt | head -n3

gunzip seq/*.gz


for i in `tail -n+2 result/metadata.txt | cut -f 1`;do
  usearch11 -fastq_mergepairs seq/${i}.R1.fq -reverse seq/${i}.R2.fq \
  -fastqout temp/${i}.merged.fq -relabel ${i}. 
done

cat temp/*.merged.fq > temp/all.fq
ls -l temp/all.fq

rm temp/*.merged.fq

gzip seq/*


usearch11 -fastx_truncate temp/all.fq \
  -stripleft 10 -stripright 10 \
  -fastqout temp/stripped.fq

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

sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa
head -n 2 temp/otus.fa

sed -i 's/\r//g' temp/otus.fa

mkdir -p result/raw
# cp temp/otus.fa result/raw/otus.fa
usearch11 -uchime2_ref temp/otus.fa -db ${db}/usearch/UNITE_modified_new2021.fasta \
     -chimeras temp/otus_chimeras.fa -strand plus -mode high_confidence -threads 30
#00:07 170Mb   100.0% Chimeras 45/2920 (1.5%), in db 248 (8.5%), not matched 2627 (90.0%)

cat temp/otus.fa temp/otus_chimeras.fa| grep '>' | sort | uniq -u \
  | sed 's/>//' > temp/non_chimeras.id

usearch11 -fastx_getseqs temp/otus.fa -labels temp/non_chimeras.id \
    -fastaout result/raw/otus.fa

 sed -i 's/\r//g' result/raw/otus.fa 

nohup usearch11 -otutab temp/filtered.fa -otus result/raw/otus.fa \
   -otutabout result/raw/otutab.txt -threads 20 &


nohup usearch11 -sintax result/raw/otus.fa -db ${db}/usearch/SILVA_138/SILVA_16S_V138.fasta \
-strand both -tabbedout result/raw/otus.sintax -sintax_cutoff 0.8 &

mv result/otutab.txt result/raw

#/mnt/data/hanyunpeng
#blastn -max_target_seqs 10 -db /data/db/NCBI_nt/nt -out aspergillus_likeout -query aspergillus_like.fas -num_threads 30 outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids salltitles"
#blast_formatter -archive "othersout" -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids salltitles" > othersout1
#blastn -max_target_seqs 10 -db /data/db/NCBI_nt/nt -out others2out -query others2.fas -num_threads 30 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids salltitles" > othersout
blastn -max_target_seqs 10 -db /data/db/NCBI_nt/nt -out penicillium5out -query penicillium5.fas -num_threads 30 -outfmt "7 qseqid pident salltitles"

nohup blastn -max_target_seqs 5 -db /data/db/NCBI_nt/nt -out Wuwei_Fungi_NCBI.blast \
-query otus.fa -num_threads 20 \
-outfmt "6 qseqid qlen qstart qend sseqid stitle slen sstart send qcovs bitscore score evalue pident" &

grep -v '^#' penicillium5out | awk '{print$1,$2,$3,$4}' > penicilliumout5 && rm -f penicillium5out

wc -l result/otutab.txt

db=../../Public_Scripts

Rscript ${db}/script/otutab_filter_nonFungi.R -h
Rscript ${db}/script/otutab_filter_nonFungi.R \
      --input result/raw/otutab.txt \
      --taxonomy result/raw/otus.sintax \
      --output result/otutab.txt\
      --stat result/raw/otutab_nonFunc.stat \
      --discard result/raw/otus.sintax.discard
      
wc -l result/otutab.txt
usearch11 -otutab_trim result/otutab.txt -min_otu_size 16 -output result/otutab1.txt
mv result/otutab1.txt result/otutab.txt

cut -f 1 result/otutab.txt | sed '1 s/#OTU ID/OTUID/' > result/otu.id
wc -l result/otu.id

usearch11 -fastx_getseqs result/raw/otus.fa -labels result/otu.id \
-fastaout result/otus.fa

awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}'\
    result/raw/otus.sintax result/otutab.id \
    > result/otus.sintax

sed -i 's/\t$/\td:Unassigned/' result/otus.sintax


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

usearch11 -otutab_stats result/otutab.txt \
  -output result/otutab.stat
cat result/otutab.stat

db=../../Public_Scripts
mkdir -p result/alpha
Rscript ${db}/script/otutab_rare.R --input result/otutab.txt \
  --depth 40907 --seed 1 \
  --normalize result/otutab_rare.txt \
  --output result/alpha/vegan.txt

usearch11 -otutab_stats result/otutab_rare.txt \
      -output result/otutab_rare.stat
cat result/otutab_rare.stat  

mkdir -p result/tax
for i in p c o f g ;do
  usearch11 -sintax_summary result/otus.sintax \
  -otutabin result/otutab_rare.txt -rank ${i} \
  -output result/tax/sum_${i}.txt
done
sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' result/tax/sum_*.txt

ls -sh result/tax/sum_*.txt
head -n3 result/tax/sum_g.txt

cat result/taxonomy2.txt | tr -s ";" "\t" > result/rep_seqs_tax.txt
head result/rep_seqs_tax.txt


#awk '{print $NF}' result/tax/test.txt
awk 'BEGIN{IFS='\t'}{$NF="";print $0}' result/tax/sum_g.txt | head -n 12 > result/tax/abundant_genus.txt

sed -i 's/[ ]*$//g' result/tax/abundant_genus.txt

#sed -i '2d' result/tax/abundant_genus.txt
cat result/tax/abundant_genus.txt | tr -s ' ' '\t' > result/tax/abundant_genus1.txt
cat result/tax/abundant_genus1.txt
mv result/tax/abundant_genus1.txt result/tax/abundant_genus.txt
cat result/tax/abundant_genus.txt

cat result/tax/final_g_group.txt | awk -F '\t' '{$2=null;print $0}' | tr -s ' ' '\t' > result/tax/final_g_group1.txt
mv result/tax/final_g_group1.txt result/tax/final_g_group.txt

mkdir -p Intercity/result
grep "Intercity" result/metadata.txt | cut -f 1 > Intercity/result/SampleID.txt

usearch11 -otutab_sample_subset result/otutab.txt \
-labels Intercity/result/SampleID.txt -output Intercity/result/otutab1.txt

usearch11 -otutab_trim Intercity/result/otutab1.txt -min_otu_size 2 -output Intercity/result/otutab.txt

usearch11 -otutab_stats Intercity/result/otutab.txt \
   -output Intercity/result/otutab.stat
cat Intercity/result/otutab.stat

cut -f 1 Intercity/result/otutab.txt | sed '1 s/#OTU ID/OTUID/' > Intercity/result/otu.id
wc -l Intercity/result/otu.id

usearch11 -fastx_getseqs result/raw/otus.fa -labels Intercity/result/otu.id \
-fastaout Intercity/result/otus.fa
mkdir -p Intercity/result/raw/

awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/otus.sintax \
Intercity/result/otu.id > Intercity/result/raw/otus.sintax
head Intercity/result/raw/otus.sintax
awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/taxonomy.txt \
Intercity/result/otu.id > Intercity/result/taxonomy.txt
cp Intercity/result/raw/otus.sintax Intercity/result/otus.sintax


usearch11 -alpha_div result/otutab_rare.txt \
  -output result/alpha/alpha.txt

#Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)-method fast  https://drive5.com/usearch/manual/cmd_otutab_subsample.html
usearch10 -alpha_div_rare result/otutab.txt \
  -output result/alpha/alpha_rare.txt -method without_replacement


usearch11 -calc_distmx result/otus.fa -tabbedout result/distmx.txt -maxdist 0.8 -termdist 0.9
  
usearch11 -otutab_core result/otutab.txt -distmxin result/distmx.txt \
  -sintaxin result/otus.sintax -tabbedout result/core.txt

Rscript ${db}/script/otu_mean.R --input result/otutab_rare.txt \
   --design result/metadata.txt \
   --group Seasonal --thre 0 \
   --output result/otutab_mean_Seasonal.txt
head -n3 result/otutab_mean_Seasonal.txt

awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
    else {for(i=2;i<=NF;i++) if($i>0.0001) print $1, a[i];}}' \
    result/otutab_mean_Seasonal.txt > result/alpha/otutab_mean_Seasonal.txt
head result/alpha/otutab_mean_Seasonal.txt


Rscript ${db}/script/alpha_rare_curve.R \
  --input result/alpha/alpha_rare.txt --design result/metadata.txt \
  --group Station --output result/alpha/ \
  --width 180 --height 120

Rscript ${db}/script/alpha_boxplot.R -h

Rscript ${db}/script/alpha_boxplot.R --alpha_index richness \
  --input result/alpha/vegan.txt --design result/metadata.txt \
  --group POSITION --output result/alpha/ \
  --width 90 --height 60

for i in `head -n1 result/alpha/vegan.txt|cut -f 2-`;do
  Rscript ${db}/script/alpha_boxplot.R --alpha_index ${i} \
     --input result/alpha/vegan.txt --design result/metadata.txt \
     --group POSITION --output result/alpha/ \
    --width 180 --height 120
done

bash ${db}/script/sp_vennDiagram.sh \
  -f result/alpha/otu_group_exist_Station_Type.txt \
  -a Intercity  -b Intercity -c Intercity \
  -w 3 -u 3 \
  -p Station_Venn


bash ${db}/script/sp_pheatmap.sh \
  -f result/beta/bray_curtis.txt \
  -H 'TRUE' -u 16 -v 12

cut -f 1,3 result/metadata.txt > temp/group.txt

bash ${db}/script/sp_pheatmap.sh \
  -f result/beta/bray_curtis.txt \
  -H 'TRUE' -u 16 -v 12 \
  -P temp/group.txt -Q temp/group.txt

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

i=o
Rscript ${db}/script/tax_circlize.R \
  --input result/tax/sum_${i}.txt --design result/metadata.txt \
  --group group --legend 10

mv circlize.pdf result/tax/sum_${i}.circlize.pdf
mv circlize_legend.pdf result/tax/sum_${i}.circlize_legend.pdf

compare="Warm-Cold"
Rscript ${db}/script/compare.R \
  --input result/otutab_rare.txt --design result/metadata.txt \
  --group Seasonal1 --compare ${compare} --threshold 0.001 \
  --method edgeR --pvalue 0.05 --fdr 0.2 \
  --output result/compare/

Rscript ${db}/script/compare_volcano.R \
   --input result/compare/${compare}.txt \
   --output result/compare/${compare}.volcano.pdf \
   --width 120 --height 120


bash ${db}/script/compare_heatmap.sh -i result/compare/${compare}.txt -l 7 \
       -d result/metadata.txt -A Seasonal1 \
       -t result/taxonomy.txt \
       -w 10 -h 14 -s 8 \
       -o result/compare/${compare}.txt

bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
   -t result/taxonomy.txt \
   -p result/tax/sum_p.txt \
   -w 180 -v 60 -s 7 -l 10 \
   -o result/compare/${compare}

bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
       -t result/taxonomy.txt \
       -p result/tax/sum_c.txt \
       -w 183 -v 80 -s 7 -l 10 -L Class \
       -o result/compare/${compare}.txt

bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
       -t result/taxonomy.txt \
       -p result/tax/sum_c.txt \
       -w 183 -v 149 -s 7 -l 10 -L Class \
       -o result/compare/${compare}.manhattan.c.legend.pdf

  Rscript ${db}/script/alpha_boxplot.R --alpha_index Lactobacillus \
    --input result/tax/sum_g.txt --design result/metadata.txt \
    --transpose TRUE --scale TRUE \
    --width 89 --height 59 \
    --group group --output result/tax/

    Rscript ${db}/script/alpha_boxplot.R --alpha_index Filobasidium \
      --input result/tax/sum_g.txt --design result/metadata.txt \
      --transpose TRUE \
      --width 89 --height 59 \
      --group group1 --output result/compare/
Rscript ${db}/script/format2stamp.R -h
mkdir -p result/stamp
Rscript ${db}/script/format2stamp.R --input result/otutab_rare.txt \
   --taxonomy result/taxonomy.txt --threshold 0.001 \
   --output result/tax

Rscript ${db}/script/format2lefse.R -h
Rscript ${db}/script/format2lefse.R --input result/otutab.txt \
  --taxonomy result/taxonomy.txt --design result/metadata.txt \
  --group StationType --threshold 0.0001 \
  --output result/LEfSe_Station_Type

mkdir -p result/compare/Ternary_plot
cut -f 2- otutab_high.mean > temp

sed '1i OTUId \t Silva_taxonomy' result/taxonomy2.txt >result/taxonomy3.txt

paste result/taxonomy3.txt result/otutab_rare.txt > result/otu_table.txt
head -n 3 result/otu_table.txt

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

mv circlize.pdf result/sum_KEGG.Pathway.circlize.pdf
mv circlize_legend.pdf result/sum_KEGG.Pathway.circlize_legend.pdf

mv result/otu_table.txt result/Ternary_plot/otu_table.txt



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

Rscript ${db}/script/otu_mean.R --input result/Pathogen_otutab_rare.txt \
   --design result/metadata.txt \
   --group Station --thre 0 \
   --output result/otutab_mean_Station.txt

awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
    else {for(i=2;i<=NF;i++) if($i>0.0001) print $1, a[i];}}' \
    result/otutab_mean_Seasonal.txt > result/alpha/otutab_mean_Seasonal.txt

cd /mnt/sdd/ZhouXin/Tomato_project/DataProcess4Publication/Field_Bacteria_16S/Bulk
mkdir -p result/OTU_tree

nohup muscle -in result/otus.fa -out result/OTU_tree/otus_aligned.fas &

trimal -in result/OTU_tree/otus_aligned.fas -out result/OTU_tree/otus_aligned_trimed.fa -gt 0.95

nohup iqtree -s result/OTU_tree/otus_aligned_trimed.fa \
        -bb 5000 -redo -alrt 5000 -nt AUTO -T 40  \
        -pre result/OTU_tree/otus &


cut -f 1 result/growthForm_sintax.txt | sed '1 s/#OTU ID/OTUID/' > result/growthForm.id

usearch11 -otutab_otu_subset result/otutab.txt -labels result/growthForm.id \
-output result/otutab_growthForm.txt

usearch11 -sintax_summary result/growthForm_sintax.txt \
  -otutabin result/otutab_growthForm.txt -rank g \
  -output result/sum_growthForm.txt

cut -f 1 result/guild_sintax.txt | sed '1 s/#OTU ID/OTUID/' > result/guild.id

usearch11 -otutab_otu_subset result/otutab.txt -labels result/guild.id \
-output result/otutab_guild.txt

usearch11 -sintax_summary result/guild_sintax.txt \
  -otutabin result/otutab_guild.txt -rank g \
  -output result/sum_guild.txt


cut -f 1 result/trait_sintax.txt | sed '1 s/#OTU ID/OTUID/' > result/trait.id

usearch11 -otutab_otu_subset result/otutab.txt -labels result/trait.id \
-output result/otutab_trait.txt

usearch11 -sintax_summary result/trait_sintax.txt \
  -otutabin result/otutab_trait.txt -rank g \
  -output result/sum_trait.txt

Rscript ${db}/script/format2lefse.R -h
Rscript ${db}/script/format2lefse.R --input result/otutab_guild.txt \
  --taxonomy Lefse/guild_taxonomy.txt --design result/metadata.txt \
  --group Station_Type --threshold 0.0001 \
  --output Lefse/Station_Type

  Rscript ${db}/script/alpha_boxplot.R --alpha_index SoftRot \
    --input result/sum_trait.txt --design result/metadata.txt \
    --transpose TRUE --scale TRUE \
    --width 89 --height 59 \
    --group group --output result/  

awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}'\
    ../otus.sintax ./Core_OTUID.txt \
    > ./otus.sintax
sed -i 's/\t$/\td:Unassigned/' ./otus.sintax


cut -f 1,4 ./otus.sintax \
  |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' \
   > ./taxonomy2.txt
head -n13 ./taxonomy2.txt

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

mkdir -p ./tax
for i in p c o f g ;do
  usearch11 -sintax_summary ./otus.sintax \
  -otutabin ./core_otutab.txt -rank ${i} \
  -output ./tax/sum_${i}.txt
done
sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' ./tax/sum_*.txt
ls -sh ./tax/sum_*.txt
cat ./taxonomy2.txt | tr -s ";" "\t" > ./rep_seqs_tax.txt
head ./rep_seqs_tax.txt

awk 'BEGIN{IFS='\t'}{$NF="";print $0}'  ./tax/sum_g.txt | head -n 15 > ./tax/abundant_genus.txt

sed -i 's/[ ]*$//g' ./tax/abundant_genus.txt

#sed -i '2d' result/tax/abundant_genus.txt
cat ./tax/abundant_genus.txt | tr -s ' ' '\t' > ./tax/abundant_genus1.txt
cat ./tax/abundant_genus1.txt
mv ./tax/abundant_genus1.txt ./tax/abundant_genus.txt
cat ./tax/abundant_genus.txt

mkdir -p result/Core_OTU_Tree
cd result/Core_OTU_Tree

usearch11 -otutab_otu_subset ../raw/otutab.txt -labels ./Core_OTUID.txt -output ./result/core_otutab.txt

usearch11 -fastx_getseqs ../otus.fa -labels ./Core_OTUID.txt \
    -fastaout ./Core_otus.fa
head -n 2 ./Core_otus.fa
grep ">" ./Core_otus.fa |wc -l
awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../taxonomy.txt \
./Core_OTUID.txt > ./Core_otutab_high.tax

awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../otutab_mean_Station_Type.txt ./Core_OTUID.txt \
    | sed 's/#OTU ID/OTUID/' > ./otutab_high_Station.mean


awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../otutab_mean_Seasonal.txt ./Core_OTUID.txt \
    | sed 's/#OTU ID/OTUID/' > ./otutab_mean_Seasonal.mean
#sed -i '2d' ./otutab_high.mean
head -n3 ./otutab_high_Station.mean

cut -f 2- ./otutab_high_Station.mean > temp1
cut -f 2- ./otutab_mean_Seasonal.mean > temp2
paste ./Core_otutab_high.tax temp1 temp2> ./annotation1.txt
head -n 3 ./annotation1.txt

cut -f 2- ./otutab_high_Station.mean > temp1
paste ./annotation1.txt temp1 > ./annotation2.txt

time muscle -in Core_otus.fa -out otus_aligned.fas

trimal -in otus_aligned.fas -out otus_aligned_trimed.fa -gt 0.95

mkdir -p iqtree
iqtree2 -s  otus_aligned_trimed.fa \
        -bb 1000 -redo -alrt 1000 -nt AUTO \
        -pre iqtree/Core_otus


Rscript ${db}/script/table2itol.R -a -c double -D plan1 -i OTUID -l Family -t %s -w 0.5 annotation.txt

Rscript ${db}/script/table2itol.R -a -d -c none -D plan2 -b Phylum -i OTUID -l Family -t %s -w 0.5 annotation.txt

Rscript ${db}/script/table2itol.R -c keep -D plan3 -i OTUID -t %s otutab_mean_Seasonal.mean

Rscript ${db}/script/table2itol.R -a -c factor -D plan4 -i OTUID -l Genus -t %s -w 0 annotation.txt

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

sed '1i OTUId \t Silva_taxonomy' result/taxonomy2.txt >result/taxonomy3.txt

paste result/taxonomy3.txt result/otutab_rare.txt > result/otu_table.txt
head -n 3 result/otu_table.txt

cat result/otu_table.txt | awk -F '\t' '{$3=null;print $0}' | tr -s ' ' '\t' > result/otu_table1.txt

sed 's/k__/D_0__/g' result/otu_table1.txt > result/otu_table.txt
sed 's/p__/D_1__/g' result/otu_table.txt > result/otu_table1.txt
sed 's/c__/D_2__/g' result/otu_table1.txt > result/otu_table.txt
sed 's/o__/D_3__/g' result/otu_table.txt > result/otu_table1.txt
sed 's/f__/D_4__/g' result/otu_table1.txt > result/otu_table.txt
sed 's/g__/D_5__/g' result/otu_table.txt > result/otu_table1.txt
sed 's/s__/D_6__/g' result/otu_table1.txt > result/otu_table.txt

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

awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/metadata.txt \
result/AutumnID.txt > result/metadata_Autumn.txt
cat temp1 result/metadata_Autumn.txt > temp_combined && mv temp_combined result/metadata_Autumn1.txt
mv result/metadata_Autumn1.txt result/metadata_Autumn.txt

grep "Winter" result/metadata.txt | cut -f 1 > result/SampleID.txt
usearch11 -otutab_sample_subset result/Pathogen_otutab_rare.txt \
-labels result/SampleID.txt -output result/Pathogen_otutab_Winter.txt

Rscript ${db}/script/otu_mean.R --input result/Pathogen_otutab_Winter.txt \
   --design result/metadata_Winter.txt \
   --group Station --thre 0 \
   --output result/otutab_mean_Station_Winter.txt
head -n3 result/otutab_mean_Station_Winter.txt
awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
    else {for(i=2;i<=NF;i++) if($i>0.0001) print $1, a[i];}}' \
    result/otutab_mean_Station.txt > result/alpha/otu_group_exist_Station.txt
head result/alpha/otu_group_exist_Station.txt


awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
    else {for(i=2;i<=NF;i++) if($i>0.0001) print $1, a[i];}}' \
    result/otutab_mean_Station_Winter.txt > result/alpha/otu_group_exist_Station_Winter.txt
head result/alpha/otu_group_exist_Station_Winter.txt

awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
    else {for(i=2;i<=NF;i++) if($i>0.0001) print $1, a[i];}}' \
    result/otutab_mean_Seasonal.txt > result/alpha/otu_group_exist_Seasonal.txt
head result/alpha/otu_group_exist_Seasonal.txt

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

  

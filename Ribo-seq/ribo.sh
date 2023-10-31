### Aligning and counting
bowtie2 -p 8 -x ~/reference/dmel/mel_rRNA -U Scr-1_R1.fastq.gz --un Scr-1_unmaped_rRNA.fasta -S Scr-1_rRNA.sam
bowtie2 -p 8 -x ~/reference/tRNA/Drosophila/tRNA_modified -U Scr-1_unmaped_rRNA.fasta --un Scr-1_unmaped_tRNA.fasta -S Scr-1_tRNA.sam
bowtie2 -p 8 -x ~/reference/dmel/dmel-longest-transcript -U Scr-1_unmaped_tRNA.fasta -S Scr-1_trans.sam
cat *_trans.sam|grep ^A|awk '{if($6!="*")print $0}'|sed 's/XG:i:1/XG:i:0/g'|sed 's/XG:i:2/XG:i:0/g'|grep "XG:i:0"|awk '{print $3"\t"$4"\t"$6}'|sort |sed 's/D/M/g'|sed 's/.I//g' > aa
cat aa|cut -f 3|sed 's/M/\t/g'|awk '{sum=0};{for (i=1;i<=NF;i++) {sum+=$i};print sum}'> cc
paste aa cc|awk '{print $1"\t"$2"\t"$4}' > Scr-1_trans.txt
join *_trans.txt /storage/SJL/reference/dmel/length_max.txt|awk '{if($2>($5+15)&& $2<($8-$7-15))print $1}'|uniq -c|awk '{print $2"\t"$1}' > Scr-1_counts.txt
###Scr-1_trans.txt, First column: FlyBase Geneid; Second column: start position; third column: length; 
###length_max.txt, First column: FlyBase Geneid; Second column: the id of transcript with longest CDS; Third column: the length of 5UTR; Fourth column: the length of CDS; Fifth column: the length of 3UTR; Sixth column: the length of transcript;  

### Determination of A site
join ~/ribo/20220607/Scr-1/Scr-1_trans.txt /storage/SJL/reference/dmel/length_max.txt|awk '{if($2>($5+15)&& $2<($8-$7-15))print $0}'|awk '{if($3>=30 && $3<=32)print $0}'|awk '{print $1"\t"$2-$5"\t"$3}'|awk '{print $1"\t"$2"\t"$2 % 3}'|awk '{if($3=="1") print $1"\t"$2+15"\t"$3}'|awk '{if($2>=0)print $1"\t"$2-1"\t"$2+2}' > Scr-1.bed
bedtools getfasta -fi /storage/SJL/ribo/20220607/Scr-1/test.fasta -bed Scr-1.bed -fo Scr-1.fasta
cat Scr-1.bed |awk '{print $1"_"$2"-"$3}'|sort -u >aa
cat Scr-1.fasta |awk '/^>/&&NR>1{print ""};{printf "%s",/^>/ ? $0" ":$0}'|sed 's/>//g'|sed 's/:/_/g'|sort -k1,1 >bb
join aa bb| sed 's/_/\t/g'|sed 's/-/\t/g'|awk '{print $4}'|sort |uniq -c|awk '{print $2"\t"$1}' > Scr-1_A.txt

### Aligning
conda activate mimseq2
cd /home/Application/anaconda3/envs/mimseq/tRNA
../../mimseq2/bin/mimseq --species Dmel --cluster --cluster-id 0.95 --threads 12 --min-cov 2000 --max-mismatches 0.1 --control-condition copy -n tRNA --out-dir test --max-multi 4 --remap --remap-mismatches 0.075 SampleData.txt

samtools view /home/Application/anaconda3/envs/mimseq/tsRNA/test4/L3-ts2_S1_L003_R1_001_trimmed.fq.gz.unpaired_uniq.bam|awk '{print $3"\t"$4"\t"$6}'|sort |sed 's/D/M/g'|sed 's/.I//g'|sed 's/.S//g' > aa
cat aa|cut -f 3|sed 's/M/\t/g'|awk '{sum=0};{for (i=1;i<=NF;i++) {sum+=$i};print sum}'> cc
paste aa cc|awk '{print $1"\t"$2"\t"$4}'|sed 's/Mrosophila_melanogaster_//g'|sed 's/Escherichia_coli_str_K_12_substr_MG1655_//g'|sed 's/tRNA-//g' |grep -v "mito"|sed 's/-/\t/g'|awk '{print $1"-"$2"\t"$5"\t"$6}'|sed 's/Rna-ATA/2SrRNA/g'|sort |sed 's/tRNAeC/SeC/g'|sed 's/tRNAer/Ser/g'|sed 's/tRNAle/Ile/g' > L3-ts1.txt

###classification
join P2-ts1.txt /storage/SJL/tsRNA/20220623/tRNA_anti.txt|awk '{if($2==1&&($2+$3-1)==$4) print $0"\t""5-tRH"; else if($2==1&&($2+$3-1)==($4+1)) print $0"\t""5-tRH";if($2==1&&($2+$3-1)==($4+2)) print $0"\t""5-tRH";else if($2==$4 && ($2+$3-1)>= $5) print $0"\t""3-tRH";else if($2==($4+1) && ($2+$3-1)>= $5) print $0"\t""3-tRH";else if ($2==($4+2) && ($2+$3-1)>= $5) print $0"\t""3-tRH";else if($2>=1&&($2+$3-1)<=($4+2)) print $0"\t""5-tRF";else if($2<=$4&&($2+$3-1)>=$4) print $0"\t""Inter-tRF";else if($2>$4) print $0"\t""3-tRF"}'|awk '{print $1"\t"$6}'|sort -k1,1|uniq -c |awk '{print $2"_"$3"\t"$1}' > L3-ts1_class.txt

### Caculating the highly expressed 7-mer tsRNA seeds
cat /storage/SJL/tsRNA/2DP/P2-1.txt| grep -v 'mt'|awk '{print $1"\t"($2-1)"\t"($2+$3-2)}' > P2-1_nc.bed
bedtools getfasta -fi /storage/SJL/reference/tRNA/Drosophila/tRNA_modified.fasta -bed P2_nc.bed -fo P2-1_nc.fasta
/storage/SJL/software/jellyfish count -m 7 -s 376929793 -t 10 P2-1_nc.fasta
/storage/SJL/software/jellyfish dump -c -t mer_counts.jf >P2-1_kmer.tsv

### After mean values of seeds expression were calculated, we used TargetScan to predict the corresponding target sites of highly expressed 7-mer tsRNA seed.
### seed_top.txt, First column: seed id; Second column: seed sequence; Third sequence: species id.
### For example: seed_1 GGUUCGA 7227
perl /storage/SJL/software/targetscan_70.pl seed_top.txt /storage/SJL/targetpredict/phyloP/fly27species_CDS.tab targetsan_CDSresults.txt
perl /storage/SJL/software/targetscan_70.pl seed_top.txt /storage/SJL/targetpredict/phyloP/fly27species_3UTR.tab targetsan_3UTRresults.txt
perl /storage/SJL/software/targetscan_70.pl seed_top.txt /storage/SJL/targetpredict/phyloP/fly27species_5UTR.tab targetsan_5UTRresults.txt

### FBtr_type.txt, First column: transcript id_seed position in mRNA.
### For example: FBtr0006151_5UTR
### input.txt, First column: transcript id_seed position in mRNA; Second column: Geneid; Third column: target position; Fourth column: target seed type; Fifth column: seed id; sixth column: FDR adjusted pvalue calcuated by phyloP; seventh column: sequence of seeds. 
### For example: FBtr0006151_5UTR       FBgn0000056     253-260 7mer-m8 seed-8  6.793841e-05    GAACCCG
for i in `cat FBtr_type.txt`
do
cat input.txt|grep $i >a
cat input.txt|grep $i|awk '{print $3}'|sed 's/-/\t/g'|awk 'NR==1{print $1;tmp=$1}NR>1{print $1-tmp;tmp=$1}'>b
paste a b|awk '{if(9!=1&&$9!=(-1)&&$9!=0)print $1"\t"$2}'|uniq -c|sed 's/_/\t/g'|awk '{print $2"\t"$4"\t"$3"\t"$1}' >>FBtr_site_num.txt
done

### FBtr_site_num.txt, First column: transcript id; Second column: Geneid; Third column: seed position in mRNA; Fourth column: site number.
### For example: FBtr0005088    FBgn0260439     3UTR    1
cat FBtr_site_num.txt|awk '{print $1}'|sort -u >test
for i in `cat test`
do
cat FBtr_site_num.txt|grep $i|awk '{sum=sum+$4};END{print $1"\t"$2"\t"sum}' >> c
done
cat c|sort -k1,1 >FBtr_num.txt

### FBtr_density.txt, First column: Geneid; Second column: transcript id; Third column: site density in mRNA. Site density in mRNA calculated by dividing the number of target sites by the transcript length.
### For example: FBgn0260439    FBtr0005088     0.00580672
### test only has one column containing Geneid.
for i in `cat test`
do
cat FBtr_density.txt|grep $i|awk 'BEGIN{max = 0}{if($3 > max) max = $3}END{print $1"\t"$2"\t"max}' >> FBtr_max_den.txt
done

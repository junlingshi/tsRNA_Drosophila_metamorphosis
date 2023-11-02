miranda /home/junling/targetpredict/5Asp.fasta /storage/SJL/reference/dmel-all-three_prime_UTR-r6.36.fasta -out 3UTR.txt
cat 3UTR.txt |grep '^>>' |sed 's/>>//g' |awk '{print $2"\t"$1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' |sort -k1,1 |awk 'BEGIN{print "Seq1""\t""Seq2""\t""Tot Score""\t""Tot Energy""\t""Max Score""\t""Max Energy""\t""Strand""\t""Len1""\t""Len2""\t""Positions"}{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > 3UTR_sorted.txt

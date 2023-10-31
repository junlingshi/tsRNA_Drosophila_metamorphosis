###FBtr_site.txt: transcriptid_startposition_endposition;
###CDS_input.txt: First column: transcriptid_startposition_endposition; Second column: species id ; Third column: sequence
for line in `cat FBtr_site.txt`
do
cat CDS_input.txt|grep $line|sort -k1,1|awk '{print ">"$2"\n"$3}'  > fail/$line".fasta"
phyloP /storage/SJL/targetpredict/phyloP/dm6.phyloP27way_modified.mod fail/$line".fasta" > fail/$line"_output.txt"
cat fail/$line"_output.txt"|grep "conservation"|sed 's/:/\t/g'|awk '{print $4}'>> pvalue
done
##FDR correction of pvalue in R language

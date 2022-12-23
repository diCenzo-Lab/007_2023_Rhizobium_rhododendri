# Make blast databases
makeblastdb -in Rhizobium_rhododendri_rho-6.2.gbk.faa -out rho62 -title rho62 -dbtype 'prot'
makeblastdb -in Rhizobium_tumorigenes_1078.gbk.faa -out 1078 -title 1078 -dbtype 'prot'
makeblastdb -in Rhizobium_tumorigenes_932.gbk.faa -out 932 -title 932 -dbtype 'prot'

# Run blast analyses
blastp -query Rhizobium_rhododendri_rho-6.2.gbk.faa -db 1078 -out rho62_1078.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send bitscore evalue sstrand' -culling_limit 1 -max_target_seqs 1 -max_hsps 1 -num_threads 28
blastp -query Rhizobium_tumorigenes_1078.gbk.faa -db rho62 -out 1078_rho62.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send bitscore evalue sstrand' -culling_limit 1 -max_target_seqs 1 -max_hsps 1 -num_threads 28
blastp -query Rhizobium_tumorigenes_1078.gbk.faa -db 932 -out 1078_932.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send bitscore evalue sstrand' -culling_limit 1 -max_target_seqs 1 -max_hsps 1 -num_threads 28
blastp -query Rhizobium_tumorigenes_932.gbk.faa -db 1078 -out 932_1078.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send bitscore evalue sstrand' -culling_limit 1 -max_target_seqs 1 -max_hsps 1 -num_threads 28

# Identify blast-BBHs
matlab -nodisplay -nosplash -nodesktop -r "find_best_hits"

# Make circos plot for rho-6.2 and 1078
cut -f1,1 -d',' orthologs.rho62.1078.txt | sed 's/ID\:Prokka\:/>ID\:Prokka\:/' | grep -v 'output' > temp1.txt
cut -f2,2 -d',' orthologs.rho62.1078.txt | sed 's/ID\:Prokka\:/>ID\:Prokka\:/' | grep -v 'output' > temp2.txt
grep -f 'temp1.txt' Rhizobium_rhododendri_rho-6.2.gbk.faa | cut -f1,7 -d'|' | sed 's/>ID\:Prokka\://' | cut -f1,2 -d':' | sed 's/\:/\ /g' | sed 's/\-/\ /g' | sed 's/|//g' | sed 's/\ /\t/g' > temp3.txt
grep -f 'temp2.txt' Rhizobium_tumorigenes_1078.gbk.faa | cut -f1,7 -d'|' | sed 's/>ID\:Prokka\://' | cut -f1,2 -d':' | sed 's/\:/\ /g' | sed 's/\-/\ /g' | sed 's/|//g' | sed 's/\ /\t/g' > temp4.txt
sed -i 's/>ID\:Prokka\://' temp1.txt
sed -i 's/>ID\:Prokka\://' temp2.txt
sed -i 's/ID\:Prokka\://g' orthologs.rho62.1078.txt
matlab -nodisplay -nosplash -nodesktop -r "merge_rho62_1078"
cat temp5.txt | sed 's/AIPMFIIK_1(3709686)/hs1/' | sed 's/AIPMFIIK_2(1530638)/hs2/' | sed 's/AIPMFIIK_3(381845)/hs3/' | sed 's/AIPMFIIK_4(336962)/hs4/' | sed 's/AMNMKMEP_1(3664408)/hs10/' | sed 's/AMNMKMEP_2(834411)/hs9/' | sed 's/AMNMKMEP_3(439071)/hs8/' | sed 's/AMNMKMEP_4(432998)/hs7/' | sed 's/AMNMKMEP_5(304572)/hs6/' | sed 's/AMNMKMEP_6(302267)/hs5/' | sed 's/,/\ /g' > links.rho62.1078.txt
rm temp*
circos -conf circos.conf.rho62.1078.txt
mv circos.png circos.rho62.1078.png
mv circos.svg circos.rho62.1078.svg

# Make circos plot for 932 and 1078
cut -f1,1 -d',' orthologs.932.1078.txt | sed 's/ID\:Prokka\:/>ID\:Prokka\:/' | grep -v 'output' > temp1.txt
cut -f2,2 -d',' orthologs.932.1078.txt | sed 's/ID\:Prokka\:/>ID\:Prokka\:/' | grep -v 'output' > temp2.txt
grep -f 'temp1.txt' Rhizobium_tumorigenes_932.gbk.faa | cut -f1,7 -d'|' | sed 's/>ID\:Prokka\://' | cut -f1,2 -d':' | sed 's/\:/\ /g' | sed 's/\-/\ /g' | sed 's/|//g' | sed 's/\ /\t/g' > temp3.txt
grep -f 'temp2.txt' Rhizobium_tumorigenes_1078.gbk.faa | cut -f1,7 -d'|' | sed 's/>ID\:Prokka\://' | cut -f1,2 -d':' | sed 's/\:/\ /g' | sed 's/\-/\ /g' | sed 's/|//g' | sed 's/\ /\t/g' > temp4.txt
sed -i 's/>ID\:Prokka\://' temp1.txt
sed -i 's/>ID\:Prokka\://' temp2.txt
sed -i 's/ID\:Prokka\://g' orthologs.932.1078.txt
matlab -nodisplay -nosplash -nodesktop -r "merge_932_1078"
cat temp5.txt | sed 's/ILJPBDBJ_1(3816680)/hs1/' | sed 's/ILJPBDBJ_2(756443)/hs2/' | sed 's/ILJPBDBJ_3(430508)/hs3/' | sed 's/ILJPBDBJ_4(339618)/hs4/' | sed 's/ILJPBDBJ_5(319653)/hs5/' | sed 's/ILJPBDBJ_6(304177)/hs6/' | sed 's/AMNMKMEP_1(3664408)/hs12/' | sed 's/AMNMKMEP_2(834411)/hs11/' | sed 's/AMNMKMEP_3(439071)/hs10/' | sed 's/AMNMKMEP_4(432998)/hs9/' | sed 's/AMNMKMEP_5(304572)/hs8/' | sed 's/AMNMKMEP_6(302267)/hs7/' | sed 's/,/\ /g' > links.932.1078.txt
rm temp*
circos -conf circos.conf.932.1078.txt
mv circos.png circos.932.1078.png
mv circos.svg circos.932.1078.svg

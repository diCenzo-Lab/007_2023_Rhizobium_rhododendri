# Manually prepared the directory repliconProteomes that includes all of the proteomes of interest
# Manually prepared the directory HMM_files that includes alignments to make the desired HMMs

# Set up directories
mkdir ProteomesHMM
mkdir HMMsearch
mkdir HMMsearchParsed
mkdir HMMsearchHits
mkdir HMMscan
mkdir hmmDatabaseFiles
mkdir HMMscanParsed
mkdir HMMscanTop
mkdir OutputFiles
mkdir HMMscanProteins
mkdir Alignments
mkdir intermediaryFiles
mkdir RAxML

# Prepare input file list
ls -1 repliconProteomes/ | sed 's/\.gbk\.faa//' > temp.txt
paste temp.txt temp.txt > intermediaryFiles/refseqGenomeInformation.txt
rm temp.txt

# Prepare protein files
perl Scripts/switchNames2.pl # switch the names of the proteins in the faa files
cat ProteomesHMM/*.faa > intermediaryFiles/combined_proteomes_HMM.faa # combine the faa files into one file
perl Scripts/modifyFasta.pl intermediaryFiles/combined_proteomes_HMM.faa > intermediaryFiles/combined_proteomes_HMM_modified.faa # modify the fasta file for easy extraction

# Perform the HMMsearch screens
# Manually added hmm_list to the directory intermediaryFiles
perl Scripts/performHMMsearch.pl # a short script to repeat for all HMM files, the build, hmmsearch, parsing, and hit extraction

# Download and prepare HMM libraries
cp /datadisk1/Users/George/1_Completed_Manuscripts/2022_northern_mesorhizobia/2_Genome_Assembly_and_Analysis/7_Sym_gene_analysis/hmmDatabaseFiles/Pfam-A.hmm.gz . # get the Pfam HMM files
cp /datadisk1/Users/George/1_Completed_Manuscripts/2022_northern_mesorhizobia/2_Genome_Assembly_and_Analysis/7_Sym_gene_analysis/hmmDatabaseFiles/TIGRFAMs_15.0_HMM.LIB.gz . # get the TIGRFAM HMM files
gunzip Pfam-A.hmm.gz # unzip the Pfam files
gunzip TIGRFAMs_15.0_HMM.LIB.gz # unzip the TIGRFAM files
mv Pfam-A.hmm hmmDatabaseFiles/Pfam-A.hmm # move the Pfam files
mv TIGRFAMs_15.0_HMM.LIB hmmDatabaseFiles/TIGRFAMs_15.0_HMM.LIB # move the TIGRFAM files
hmmconvert hmmDatabaseFiles/Pfam-A.hmm > hmmDatabaseFiles/Pfam-A_converted.hmm # convert the database to the necessary format
hmmconvert hmmDatabaseFiles/TIGRFAMs_15.0_HMM.LIB > hmmDatabaseFiles/TIGRFAM_converted.hmm # convert the database to the necessary format
cat hmmDatabaseFiles/Pfam-A_converted.hmm hmmDatabaseFiles/TIGRFAM_converted.hmm > hmmDatabaseFiles/converted_combined.hmm # combined all hidden Markov models into a single file
hmmpress hmmDatabaseFiles/converted_combined.hmm # prepare files for hmmscan searches

# Perform the HMM scan screens
perl Scripts/performHMMscan.pl # a short script to repeat for all the HMM search output files, to perform hmmscan, parse, and hit extraction

# Get hits from HMMscan screens
perl Scripts/parseHMMscanHits.pl # a script to extract the names of proteins that have top hits to the desired HMMs
LC_ALL=C sort -k3,3 -g HMMscanTop/ParA_names.txt > HMMscanTop/ParA_names_sorted.txt
LC_ALL=C sort -k3,3 -g HMMscanTop/ParBc_names.txt > HMMscanTop/ParBc_names_sorted.txt
LC_ALL=C sort -k3,3 -g HMMscanTop/Rep_3_names.txt > HMMscanTop/Rep_3_names_sorted.txt
sort -k2,2 HMMscanTop/ParA_names_sorted.txt > HMMscanTop/ParA_names_unique.txt
sort -k2,2 HMMscanTop/ParBc_names_sorted.txt > HMMscanTop/ParBc_names_unique.txt
sort -k2,2 HMMscanTop/Rep_3_names_sorted.txt > HMMscanTop/Rep_3_names_unique.txt
perl Scripts/extractHMMscanHits.pl HMMscanTop/ParA_names_unique.txt > HMMscanProteins/ParA.faa # extract the ParA proteins
perl Scripts/extractHMMscanHits.pl HMMscanTop/ParBc_names_unique.txt > HMMscanProteins/ParBc.faa # extract the ParB proteins
perl Scripts/extractHMMscanHits.pl HMMscanTop/Rep_3_names_unique.txt > HMMscanProteins/Rep_3.faa # extract the Rep proteins

# Align and trim the hits
mafft --localpair --thread 28 HMMscanProteins/ParA.faa > Alignments/ParA_mafft.faa # align ParA proteins
mafft --localpair --thread 28 HMMscanProteins/ParBc.faa > Alignments/ParBc_mafft.faa # align ParB proteins
mafft --localpair --thread 28 HMMscanProteins/Rep_3.faa > Alignments/Rep_3_mafft.faa # align Rep proteins
trimal -in Alignments/ParA_mafft.faa -out Alignments/ParA_trimal.faa -fasta -automated1 # trim ParA alignment
trimal -in Alignments/ParBc_mafft.faa -out Alignments/ParBc_trimal.faa -fasta -automated1 # trim ParB alignment
trimal -in Alignments/Rep_3_mafft.faa -out Alignments/Rep_3_trimal.faa -fasta -automated1 # trim Rep alignment

# Make fasttree phylogeny
fasttree Alignments/ParA_trimal.faa > OutputFiles/ParA_fasttree.tre # phylogeny of ParA
fasttree Alignments/ParBc_trimal.faa > OutputFiles/ParBc_fasttree.tre # phylogeny of ParB
fasttree Alignments/Rep_3_trimal.faa > OutputFiles/Rep_3_fasttree.tre # phylogeny of Rep

# Make RAxML phylogeny
cd RAxML # change directory
cp ../Alignments/ParA_trimal.faa . # get the input alignment
raxmlHPC-HYBRID-AVX2 -T 10 -s ParA_trimal.faa -N 5 -n test_phylogeny -f a -p 12345 -x 12345 -m PROTGAMMAAUTO # figure out which model to use
mpiexec --map-by node -np 5 raxmlHPC-HYBRID-AVX2 -T 6 -s ParA_trimal.faa -N autoMRE -n ParA_phylogeny -f a -p 12345 -x 12345 -m PROTGAMMALGF # make ML phylogeny
cd .. # change directory
cp RAxML/RAxML_bipartitions.ParA_phylogeny OutputFiles/ParA_raxml.tre # copy final phylogeny to output folder

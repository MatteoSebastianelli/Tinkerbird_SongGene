#This script contains the commands used to BLAST the tinkerbird reference genome against the zebra finch reference genome.

#1. shorten the NCBI headers from the two genomes
sed '/^>/ s/ .*//' GCA_015220805.1_bPogPus1.pri_genomic.fna > GCA_015220805.1_bPogPus1.pri_genomic_renamed.fna #tinkerbird
sed '/^>/ s/ .*//' GCF_003957565.2_bTaeGut1.4.pri_genomic.fna > GCF_003957565.2_bTaeGut1.4.pri_genomic_renamed.fna #zebra finch

#2. generate the zebra finch database for BLAST
makeblastdb -in GCF_003957565.2_bTaeGut1.4.pri_genomic_renamed.fna -parse_seqids -dbtype nucl -out zebrafinch_db -logfile zebrafinch_makeblastdb.log
 
#3. extract the fasta sequences of the tinkerbird genome containing the significan GEMMA SNPs

#first_snps.bed
#chr25	5000000	12000000
#last_snp.bed
#chr25	17000000	18000000

bedtools getfasta -fi GCA_015220805.1_bPogPus1.pri_genomic_renamed.fna -bed first_snps.bed -fo first_snps.fasta
bedtools getfasta -fi GCA_015220805.1_bPogPus1.pri_genomic_renamed.fna -bed last_snp.bed -fo last_snp.fasta

#3. Blast the two fasta files against the zebra finch genome
blastn -num_threads 32 -task dc-megablast -outfmt 6 -db zebrafinch_db -query first_snps.fasta -out BLAST_first_snps.out 
blastn -num_threads 32 -task dc-megablast -outfmt 6 -db zebrafinch_db -query last_snp.fasta -out BLAST_last_snp.out 

#4. Remove small matches (< 200 bp)
awk '(NR>1) && ($4 > 200)' BLAST_first_snps.out > BLAST_first_snps_no_small.out 
awk '(NR>1) && ($4 > 200)' BLAST_last_snp.out > BLAST_last_snp_no_small.out 

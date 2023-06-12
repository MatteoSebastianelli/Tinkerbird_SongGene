################################################################################
					SNP CALLING PIPELINE (PAIRED-END READS, ddRAD)
################################################################################

################################################################################
						STEP 1 - PRE-PROCESSING
################################################################################

___________________________1.1_Process_RadTag_STACKS_v2.62____________________________________
## Deplexing barcodes
## Tab-delimited Barcode file specifies which sample has which barcode, A1 -> A12; B1 -> B12; etc until H1 -> H12, col 1: barcode sequence & col 2: sample ID
## output is that every pooled sample is now represented by two files: read1 and read2
# -1 specifies forward reads
# -2 specifies reverse reads 
# -b specifies barcodes file
# -i species input file type
# -e specifies the enzyme 
# -r rescues barcodes (default is likely allowing 2 mismatches)
# -c removes uncalled bases 
# -q removes low quality (default Seq Quality <10)
# -o specifies output directory
# -y species output file type 

for file in *_read1.fq.gz
do
	process_radtags \
	-1 $file \
	-2 ${file:0:-12}_read2.fq.gz \
	-b ${file:0:-12}.txt \
	-i gzfastq \
	-e sbfI \
	-r \
	-c \
	-o /demultiplexed_reads \
	-y gzfastq
done


___________________________1.2_Clone_Filter_STACKS_v2.62____________________________________
## Removes PCR duplicates
## Last place where paired reads are needed
## Specify read pairs by hand (-1 & -2) 
## output is another pair of files for each sample with the naming: prefix.1.1.fq

for sample in *.1.fq.gz
do
	clone_filter \
	-1 $sample \
	-2 ${sample:0:-8}.2.fq.gz \
	-i gzfastq \
	-o /clone_filtered \
	-y gzfastq
done

################################################################################
						STEP 2 - MAPPING/ALIGNMENT
################################################################################

________________________2a_BWA_v0.7.1_and_2b_Samtools_v1.9______________________________
##Align and map reads to reference genome:
#-t=threads 
#-M=Mark shorter split hits as secondary (for Picard compatibility)
#-p=interleaved fq file 
#-R=creates defined read groups (clearly defined read groups are needed for BQSR)

##Convert aligned/mapped sam file into bam file, sort it by coordinate, and index it:
#-@=threads 
#-b=bam output file 
#-q=skip alignments with MAPQ smaller than defined integer
#-O=bam output file

for file in *.1.1.fq.gz
do	
	bwa mem \
		-t 2 -M \
		bPogPus1.pri.cur.20200514.fasta \
		$file \
		${file:0:-10}.2.2.fq.gz | 
	samtools view \
		-@ 2 -bS -q 20 |
	samtools sort \
		-@ 2 -O bam \
		-o /aligned_sorted/${file:0:-10}_sorted.bam
done


################################################################################
						STEP 3 - Call SNPs
################################################################################

___________________________3.1_gstacks_STACKS_v2.62____________________________________
##Identify SNPs and genotype each individuals within a specified population
#-O specifies the output directory (directory must be pre-made)
#-I specifies the directory containing aligned and sorted bams
#-M specifies the tab seperated population map with list of samples and their corresponding population (sample_name	pop)
#--var-alpha specifies the alpha threshold for discovering SNPs (smaller = more stringent)
#--gt-alpha specifies the alpha threshold for calling genotypes (smaller = more stringent)
#-t specifies the number of threads

gstacks \
-O gstacks_population1 \
-I aligned_sorted \
-M population1_popmap.txt \
--var-alpha 0.01 \
--gt-alpha 0.01 \
-t 4


___________________________3.2_populations_STACKS_v2.62____________________________________
##Calculate population genetics statistics and export to vcf format
#-P specifies the directory containing the gstacks file
#-O specifies the output directory (directory must be pre-made)
#-M specifies the tab seperated population map with list of samples and their corresponding population (sample_name	pop)
#--vcf specifies the output format
#-t specifies the number of threads

populations \
-P gstacks_population1 \
-O populations_population1 \
-M population1_popmap.txt \
--vcf \
-t 4
################################################################################
					SNP CALLING PIPELINE (PAIRED-END READS)
################################################################################


################################################################################
						STEP 1 - PRE-PROCESSING
################################################################################

___________________________1.1_CutAdapt_v4.0____________________________________
#Remove adaptor sequences from paired ends and write to a single file:
#--interleaved = writes files to a single interleaved file
#-q = trims bases lower than defined quality from 3' end of both pairs before trimming adaptors
#-a = 3' adaptor read 1; -A = 3'adaptor read 2 
#-m = discards reads shorter than defined size

conda activate cutadaptenv

for file in *_1.fq.gz
do
cutadapt \
	--interleaved \
	-q 20 \
	-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
	-m 20 \
	-o ${file:0:-8}.trm.int.fq.gz \
	$file \
	${file:0:-8}_2.fq.gz
done


________________________1.2a_BWA_v0.7.1_and_1.2b_Samtools_v1.9______________________________
#Align and map reads to reference genome:
#-t=threads 
#-M=Mark shorter split hits as secondary (for Picard compatibility)
#-p=interleaved fq file 
#-R=creates defined read groups (clearly defined read groups are needed for BQSR)

#Convert aligned/mapped sam file into bam file, sort it by coordinate, and index it:
#-@=threads 
#-b=bam output file 
#-q=skip alignments with MAPQ smaller than defined integer
#-O=bam output file

for file in *.fq.gz
do	
	bwa mem \
		-t 4 -M -p \
		-R "@RG\\tID:${file:0:7}\\tPL:illumina\\tSM:${file:0:7}\\tPU:${file:25:-14}\\tLB:${file:8:-30}" \
		bPogPus1_combined_assembly.fasta \
		$file |
	samtools-1.9/samtools view \
		-@ 4 -b -q 20 |
	samtools sort \
		-@ 4 -O bam \
		-o ${file:0:7}_${file:25:-17}.mapped.sorted.bam
	samtools index ${file:0:7}_${file:25:-17}.mapped.sorted.bam
done



________________________1.3_Picard_v2.23.1_Mark_Duplicates______________________________
#Mark duplicate reads 
#-AS = assume that the data is coordinate sorted

#Two options:
#1.) Can merge same sample from multiple runs (saved as sampleID_Run.mapped.sorted.bam)
#2.) If sample only has one run, then use single sample only

for file in *_HCHJ7DSXY.mapped.sorted.bam
do
	((i=i%N)); ((i++==0)) && wait
	java -Xmx4g -jar picard.jar MarkDuplicates \
		I=$file \
		I=${file:0:-28}_HGYLYDSXY.mapped.sorted.bam \
		O=${file:0:7}_1.merged.dupsrem.bam \
		M=${file:0:7}_1.merged.dupsrem.met.txt \
		AS=true \
		REMOVE_DUPLICATES=true \
		CREATE_INDEX=true &
done

for file in *.mapped.sorted.bam
do
	((i=i%N)); ((i++==0)) && wait
	java -Xmx4g -jar picard.jar MarkDuplicates \
		I=$file \
		O=${file:0:7}_2.merged.dupsrem.bam \
		M=${file:0:7}_2.merged.dupsrem.met.txt \
		AS=true \
		REMOVE_DUPLICATES=true \
		CREATE_INDEX=true &
done



________________________1.4_Samtools_v1.9______________________________
#Merge samples that have 3 or more runs 
#(picard mark dups can only merge 2 runs at a time)

for file in *_1.merged.dupsrem.bam
do
	samtools merge \
	-O bam \
	-@ 4 \
	${file:0:-21}.merged.dupsrem.bam \ #this is the output file name
	$file \
	${file:0:-21}_2.merged.dupsrem.bam 
	samtools index ${sample:0:-21}.merged.dupsrem.bam
done



________________________1.5_Samtools_v1.9______________________________
#Get read depth of each sample

for file in *.merged.dupsrem.bam
do
	samtools depth $file
done


##########################################################################
############## GATK 4 No Longer Uses RealignerTargetCreator ##############
### According to GATK recommended best practices, not necessary to use ###
###### IndelRealigner if you're using HaplotypeCaller to call SNPS #######
##########################################################################




################################################################################
				STEP 2 - PRE-BQSR - Create "know_variants" file
################################################################################

______2.1a_GATK_v4.0.11.0____2.1b_VCFTools_v0.1.16___2.1c_htslib_v1.9___________
#Call SNPs for individual samples for each reference genome interval
#Merge each interval 
#BGZIP and Index merged sample

sample_list=($(<samples.list))
sample=${sample_list[${SLURM_ARRAY_TASK_ID}]}

intervallist=("intervals1.list"
			"intervals2.list"
			"intervals3.list"
			"intervals4.list"
			"intervals5.list"
			"intervals6.list"
			"intervals7.list"
			"intervals8.list"
			"intervals9.list"
			"intervals10.list"
			"intervals11.list"
			"intervals12.list"
			"intervals13.list"
			"intervals14.list"
			"intervals15.list"
			"intervals16.list"
			"intervals17.list"
			"intervals18.list"
			"intervals19.list"
			"intervals20.list")

for interval in ${list[@]}
do
	java -Xmx4g -jar gatk HaplotypeCaller \
	-R bPogPus1_combined_assembly.fasta \
	-I $sample \
	-L $interval \
	-ERC GVCF \
	-O ${sample:0:-10}_${interval:0:-5}.g.vcf
done

vcf-concat \
${sample:0:-10}_intervals1.g.vcf \
${sample:0:-10}_intervals2.g.vcf \
${sample:0:-10}_intervals3.g.vcf \
${sample:0:-10}_intervals4.g.vcf \
${sample:0:-10}_intervals5.g.vcf \
${sample:0:-10}_intervals6.g.vcf \
${sample:0:-10}_intervals7.g.vcf \
${sample:0:-10}_intervals8.g.vcf \
${sample:0:-10}_intervals9.g.vcf \
${sample:0:-10}_intervals10.g.vcf \
${sample:0:-10}_intervals11.g.vcf \
${sample:0:-10}_intervals12.g.vcf \
${sample:0:-10}_intervals13.g.vcf \
${sample:0:-10}_intervals14.g.vcf \
${sample:0:-10}_intervals15.g.vcf \
${sample:0:-10}_intervals16.g.vcf \
${sample:0:-10}_intervals17.g.vcf \
${sample:0:-10}_intervals18.g.vcf \
${sample:0:-10}_intervals19.g.vcf \
${sample:0:-10}_intervals20.g.vcf \
> ${sample:0:-10}_combined.g.vcf

bgzip ${sample:0:-10}_combined.g.vcf
tabix -s1 -b2 -e2 ${sample:0:-10}_combined.g.vcf.gz



____________________________2.2_GATK_v4.0.11.0___________________________________
#Combine all samples for subsequent joint genotyping
#List each sample using the "--variant" flag 
#(unfortunately, does not use a list of samples)

java -Xmx16g -jar gatk CombineGVCFs \
	-R bPogPus1_combined_assembly.fasta \
	--variant sample1_combined.g.vcf.gz \
	--variant sample2_combined.g.vcf.gz \
	--variant samlpe3_combined.g.vcf.gz \
	--variant sample4_combined.g.vcf.gz \
	--variant sample5_combined.g.vcf.gz \
	-O all_samples_combined_pre_bqsr.g.vcf.gz



____________________________2.3_GATK_v4.0.11.0___________________________________
#Joint genotype

java -Xmx4g -jar gatk GenotypeGVCFs \
		-R bPogPus1_combined_assembly.fasta \
		-V all_samples_combined_pre_bqsr.g.vcf.gz \
		-O all_samples_combined_pre_bqsr_genotyped.vcf.gz



___________________2.4a_BCFTools_v1.9___2.4b_htslib_v1.9__________________________
#Use GATK recommended hard filter parameters to create 
#files with variants which will be used for BQSR

#this has to be created using our samples as this is a non-model organism
#model organisms usually have a pre-made file of variants known to be in that species

bcftools view \
-m2 -M2 -v snps \
all_samples_combined_pre_bqsr_genotyped.vcf.gz | 
bcftools filter \
-i 'QUAL>30 && FS<60.0 && SOR<3 && MQ>40.0 && MQRankSum>-12.5 && MQRankSum<12.5 && QD>2.0 && ReadPosRankSum>-8.0 && ReadPosRankSum<8.0' |
bcftools filter \
-e 'FMT/DP<3 | FMT/GQ<20' |
bcftools filter \
-e 'AC==0 || AC==AN' \
--SnpGap 10 \
-Oz -o known_snps.vcf.gz

tabix -s1 -b2 -e2 known_snps.vcf.gz

bcftools view \
-v indels \
all_samples_combined_pre_bqsr_genotyped.vcf.gz | 
bcftools filter \
-i 'QUAL>30 && FS<200.0 && QD>2.0 && ReadPosRankSum>-20.0' \
-Oz -o known_indels.vcf.gz

tabix -s1 -b2 -e2 known_indels.vcf.gz



################################################################################
						STEP 3 - BQSR
################################################################################

____________________________3.1-3.4_GATK_v4.0.11.0______________________________
#Create bqsr tables which will be used to produce graphs 
#of the changes made in samples after applying bqsr

#Use pre-processed bam file and the "known_variants" files to create table
for bam in *.bam
do
	java -Xmx4g -jar gatk BaseRecalibrator \
	-I $bam \
	-R bPogPus1_combined_assembly.fasta \
	-O ${bam:0:-19})_recal1.table \
	-known-sites known_snps.vcf.gz \
	-known-sites known_indels.vcf.gz
done


#Use recalibration table from previous step 
for bam in *.bam
do
	java -Xmx4g -jar gatk ApplyBQSR \
	-I $bam \
	-R bPogPus1_combined_assembly.fasta \
	-O ${bam:0:-19}_recal.bam \
	--bqsr-recal-file ${bam:0:-19}_recal1.table
done


#Use newly recalibrated bam files to create another table
for bam in *_recal.bam;
do
	java -Xmx4g -jar gatk BaseRecalibrator \
	-I $bam \
	-R bPogPus1_combined_assembly.fasta \
	-O ${bam:0:-4}2.table \
	-known-sites known_snps.vcf.gz \
	-known-sites known_indels.vcf.gz
done


#Use the two tables to make graphs
for table in *1.table
do
    java -Xmx4g -jar gatk AnalyzeCovariates \
     -before $table \
     -after ${table:0:-7}2.table \
     -plots ${table:0:-7}plots.pdf
done



################################################################################
				STEP 4 - POST-BQSR - Call SNPs for final vcf
################################################################################

____________________________4.1-4.4_GATK_v4.0.11.0_______________________________
#Using the recalibrated bams
#Repeat steps 2.1 (HaplotyeCaller), 
#2.2 (CombineGVCFs), and
#2.3 (GenotypeVCFs) to produce final vcf
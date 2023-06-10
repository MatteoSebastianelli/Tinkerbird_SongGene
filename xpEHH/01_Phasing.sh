##########################################################

# Statistical phasing
# This script follows Mark Ravinet's speciation genomics tutorial (https://speciationgenomics.github.io/phasing/)

#tinkerbirds_merged_SNPs_biallelic_max-miss0.95_maf0.03_minmaxDP3.5_50.vcf:
#miss = 0.95, 
#maf = 0.03, 
#minDP = 3.5
#maxDP = 50

##########################################################

# Extract chr. 25 only for allopatric individuals (listed in the AllopatricPop.txt file)
vcftools --vcf tinkerbirds_merged_SNPs_biallelic_max-miss0.95_minmaxDP3.5_50.vcf --chr SUPER_25 --keep AllopatricPop.txt --recode --recode-INFO-all --out tinkerbirds_alloatric_chr25_SNPs_biallelic_max-miss0.95_minmaxDP3.5_50

# Set input and output variables
VCF=./tinkerbirds_alloatric_chr25_SNPs_biallelic_max-miss0.95_minmaxDP3.5_50.vcf

OUTPUT=${VCF%_u*}_phased

# Run shapeit
shapeit --input-vcf $VCF -O $OUTPUT --window 0.5 -T 4

# Convert shapeit output to vcf
shapeit -convert --input-haps ${OUTPUT} --output-vcf ${OUTPUT}.vcf

# Compress vcf file
bgzip ${OUTPUT}.vcf

# Index vcf file
bcftools index ${OUTPUT}.vcf.gz

# Make a list of the individuals in each population i.e. pop_YFTallo.txt and pop_RFTallo.txt
# First, set a new variable
VCF=./tinkerbirds_alloatric_chr25_SNPs_biallelic_max-miss0.95_minmaxDP3.5_50_phased.vcf.gz

# Index vcf file
/tigress/VONHOLDT/Simona/miniconda3/bin/bgzip $VCF

# Split each vcf in the two populations 
bcftools view -S pop_YFTallo.txt -O z -o Chr25_yft.vcf.gz $VCF

bcftools view -S pop_RFTallo.txt -O z -o Chr25_rft.vcf.gz $VCF
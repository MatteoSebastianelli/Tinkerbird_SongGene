##########################################################

# Run association test in GEMMA

#tinkerbirds_merged_SNPs_biallelic_max-miss0.95_maf0.03_minmaxDP3.5_50.vcf.gz:
#miss = 0.95, 
#maf = 0.03, 
#minDP = 3.5
#maxDP = 50

##########################################################
#Extract recorded individuals
plink --vcf tinkerbirds_merged_SNPs_biallelic_max-miss0.95_maf0.03_minmaxDP3.5_50.vcf --allow-extra-chr --allow-no-sex --keep KeepFor_songgwas_sympatry.txt --make-bed --out tinkerbirds_ioigwas_SNPs_biallelic_max-miss0.95_maf0.03_minmaxDP3.5_50

# Generate relatedness matrix
gemma-0.98.5-linux-static-AMD64 -bfile tinkerbirds_ioigwas_SNPs_biallelic_max-miss0.95_maf0.03_minmaxDP3.5_50 -gk 1 -o gemma_relatedness_ioi_Symp

# Run association test
gemma-0.98.5-linux-static-AMD64 -bfile tinkerbirds_ioigwas_SNPs_biallelic_max-miss0.95_maf0.03_minmaxDP3.5_50 -k gemma_relatedness_ioi_Symp.cXX.txt -lmm 4 -o gemma_ioi_Symp_max-miss0.95_maf0.03_minmaxDP3.5_50
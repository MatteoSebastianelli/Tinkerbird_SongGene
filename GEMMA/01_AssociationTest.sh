##########################################################

# Run association test in GEMMA

#tinkerbirds_merged_SNPs_biallelic_max-miss0.95_maf0.03_minmaxDP3.5_50:
#miss = 0.95, 
#maf = 0.03, 
#minDP = 3.5
#maxDP = 50

##########################################################

# Generate relatedness matrix
gemma-0.98.5-linux-static-AMD64 -bfile tinkerbirds_merged_SNPs_biallelic_max-miss0.95_maf0.03_minmaxDP3.5_50 -gk 1 -o gemma_relatedness_ioi_Symp

# Run association test
gemma-0.98.5-linux-static-AMD64 -bfile tinkerbirds_merged_SNPs_biallelic_max-miss0.95_maf0.03_minmaxDP3.5_50 -k gemma_relatedness_ioi_Symp.cXX.txt -lmm 4 -o gemma_ioi_Symp_max-miss0.95_maf0.03_minmaxDP3.5_50
##########################################################

# Bayesian sparse linear mixed model (BSLMM)

#tinkerbirds_merged_SNPs_biallelic_max-miss0.95_maf0.03_minmaxDP3.5_50:
#miss = 0.95, 
#maf = 0.03, 
#minDP = 3.5
#maxDP = 50

#The BSLMM will predict the phenotype for the individuals whose phenotype is NA (or -9).
#Create as many .fam file copies as the number of individuals. Each .fam file will have
#the phenotype of the individual of interest set to -9, so that the information of the other
#individuals will be used to predict its phenotype. I keep all these .fam files in a separate folder, here called fam_files

##########################################################

for ind in `cat List_87RecordedTinkerbirds.txt` # List of individuals
do

# Move fam files in main folder and rename it like the other PLINK files
cp ./fam_files/tinkerbirds_ioigwas_SNPs_biallelic_max-miss0.95_maf0.03_minmaxDP3.5_50"$ind".fam /scratch/gpfs/ms0553/matteo/gemma/bslmm/tinkerbirds_ioigwas_SNPs_biallelic_max-miss0.95_maf0.03_minmaxDP3.5_50.fam

# Compute hyper-parameters
gemma -bfile tinkerbirds_ioigwas_SNPs_biallelic_max-miss0.95_maf0.03_minmaxDP3.5_50 -bslmm 1 -k gemma_relatedness_ioi_Symp.cXX.txt -w 5000000 -s 20000000 -rpace 1000 -o "$ind"_bslmm

# Predict phenotypes using output from BSLMM
gemma -bfile tinkerbirds_ioigwas_SNPs_biallelic_max-miss0.95_maf0.03_minmaxDP3.5_50 -epm ./output/"$ind"_bslmm1.param.txt -emu ./output/"$ind"_bslmm.log.txt -ebv ./output/"$ind"_bslmm.bv.txt -k gemma_relatedness_ioi_Symp.cXX.txt -predict 1 -o "$ind"_bslmm_pred

rm tinkerbirds_ioigwas_SNPs_biallelic_max-miss0.95_maf0.03_minmaxDP3.5_50.fam

done
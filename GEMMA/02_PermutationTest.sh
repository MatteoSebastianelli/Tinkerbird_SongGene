##########################################################

# Permutation test

#tinkerbirds_ioigwas_SNPs_biallelic_max-miss0.95_maf0.03_minmaxDP3.5_50:
#miss = 0.95, 
#maf = 0.03, 
#minDP = 3.5
#maxDP = 50

##########################################################

for i in {001..1000}
do


#Randomize phenotypes:
#from the plink .fam file create a pheno.txt file with just the phenotypes and 
#a nopheno.txt file with just the first 5 columns of the .fam file

#For every loop, shuffle the phenotypes and merge them with the nopheno.txt file. Make sure to rename the resulting .fam file
# as the .bed and .bim file for the association test
shuf tinkerbids_pheno.txt | paste nopheno.fam - > tinkerbirds_ioigwas_SNPs_biallelic_max-miss0.95_maf0.03_minmaxDP3.5_50.fam

#Run gemma linear mixed model
gemma-0.98.5-linux-static-AMD64 -bfile tinkerbirds_ioigwas_SNPs_biallelic_max-miss0.95_maf0.03_minmaxDP3.5_50 -k gemma_relatedness_ioi_Symp.cXX.txt -lmm 4 -o "$i"gemma_ioi_Symp_max-miss0.95_maf0.03_minmaxDP3.5_50

#Record summary output
cd output

cut -f 13 "$i"gemma_ioi_Symp_max-miss0.95_maf0.03_minmaxDP3.5_50.assoc.txt > "$i"waldpval.txt
cat "$i"waldpval.txt | awk  'BEGIN{min=1}{if(($1)<min)  min=($1)}END {print min}' >> "$i"smallestpval.txt 
cat "$i"waldpval.txt | awk '$1<1e-6' | wc -l >> "$i"snpsbelow1e-6.txt
sort -n "$i"waldpval.txt | grep -v "nan" | grep -v "p_wald" | head -10 > "$i"_10smallestP.txt

done
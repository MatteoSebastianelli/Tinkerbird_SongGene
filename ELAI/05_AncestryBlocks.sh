##########################################################

# Estimate ancestry blocks

##########################################################

for ind in `cat /List_ToInfer.txt`
do

cd /"$ind"

# 1) Convert ancestry estimates in categories:

	R CMD BATCH /scratch/gpfs/ms0553/matteo/ELAI/MS_ELAI_ancestry_final_v3.R # ancestry score <= 0.5 = homozygous extoni
       																			# ancestry score > 0.5 & < 1.5 = heterozygous
       																			# ancestry score >= 1.5 = homozygous pusillus

    cat Chr*_combined_ave_anc_pos.txt > "$ind"_AllChr_combined_ave_anc_pos.txt # merge all chr

# 2) Run ancestryblocks_v4.pl in loop to calculate ancestry blocks

    for file in /"$ind"/Final_*.txt
    do
    	perl ./ancestryblocks_v4.pl "$file" 

    done

# Calculating genome-wide ancestry

       cat Final_*_blocks.txt > Final_"$ind"_allChr_blocks.txt 

       awk -v ind="$ind" '{print ind, $1, $2, $3, $5, $4}' FS='\t' OFS='\t' Final_"$ind"_allChr_blocks.txt > Final_"$ind"_allChr_blocks_ancestry.txt

       perl /scratch/gpfs/ms0553/matteo/ELAI/ancestrycalulator.pl < Final_"$ind"_allChr_blocks_ancestry.txt > "$ind"_WGancestry.txt

       mv *WGancestry.txt /WG_ancestry # Move output to final folder


done
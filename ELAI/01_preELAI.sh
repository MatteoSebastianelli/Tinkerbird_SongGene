##########################################################

# running PLINK on the populations output to get .geno.txt
# ad .pos.txt PLINK files

#tinkerbirds138_ELAI_ready:
#miss = 0.95, 
#maf = 0.03, 
#minDP = 3.5
#maxDP = 50

##########################################################
## make a "keep.txt" file for each reference population and for the individuals you will estimate ancestry for
## make sure you use the right flag if your system has more chr than humans
## tinkerbirds have 44 autosomes plus W and Z (which i have renamed as 45 and 46)
## Keep all associated outfiles in the direction where you will run ELAI

##Parental group 1: extoni
for chr in {01..46}
do

plink --bfile tinkerbirds138_ELAI_ready --allow-no-sex --chr-set 44 no-xy --allow-extra-chr 0 --missing-genotype 0 --chr "$chr" --recode bimbam --keep keep_for_elai_extoni.txt --out ref_extoni_for_elai_chr"$chr"

done


##Parental group 2: pusillus
for chr in {01..46}
do

plink --bfile tinkerbirds138_ELAI_ready --allow-no-sex --chr-set 44 no-xy --allow-extra-chr 0 --missing-genotype 0 --chr "$chr" --recode bimbam --keep keep_for_elai_pusillus.txt --out ref_pusillus_for_elai_chr"$chr"

done

## Group to infer ancestry
for chr in {01..46}
do

plink_linux_x86_64/plink --bfile tinkerbirds138_ELAI_ready --allow-no-sex --chr-set 44 no-xy --allow-extra-chr 0 --missing-genotype 0 --chr "$chr" --recode bimbam --keep keep_for_elai_to_infer.txt --out chrysopus_to_infer_for_elai_chr"$chr"

done




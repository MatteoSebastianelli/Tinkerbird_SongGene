##########################################################

# Combine and average ELAI runs

##########################################################

mkdir chrom_position

for chr in {01..46}
do
    chrpos_elai=/scratch/gpfs/ms0553/matteo/ELAI/elai_output
    awk '{print $3 "\t" $2}' $chrpos_elai/chrysopus_to_infer_for_elai_chr"$chr".recode.pos.txt > chrom_position/Chr"$chr"_Positions.txt # Use .pos.txt files to get chromosome positions
done


for ind in `cat List_ToInfer.txt` #This is a list of the IDs of the individuals whose ancestry needs to be inferred
do

    mkdir $ind #Create individual folders

##################################Step 1##################################
##############################Combine files###############################

	for chr in {01..46}
    do

        cat runs/chr"$chr"*"$ind".txt > $ind/Chr"$chr"."$ind"_combined.txt


##################################Step 2##################################
##############################Average files###############################
##average every column in the combined file 

        awk '{for(i=1; i<=NF; i++) {a[i]+=$i; if($i!="") b[i]++}}; END {for(i=1; i<=NF; i++) printf "%s%s", a[i]/b[i], (i==NF?ORS:OFS)}' $ind/Chr"$chr"."$ind"_combined.txt > $ind/Chr"$chr"."$ind"_combined_ave.txt


##################################Step 3##################################
################keep just one ancestry and transpose #####################

        awk '{ for (i=2;i<=NF;i+=2) print $i }' $ind/Chr"$chr"."$ind"_combined_ave.txt > $ind/Chr"$chr"."$ind"_combined_ave_anc.txt

####add chr positions
##################################Step 4##################################
###########################add chr positions##############################

        awk -v ind="$ind" '{print ind}' $ind/Chr"$chr"."$ind"_combined_ave_anc.txt > $ind/Chr"$chr"."$ind"_ID_length.txt

		paste -d "\t" chrom_position/Chr"$chr"_Positions.txt $ind/Chr"$chr"."$ind"_combined_ave_anc.txt $ind/Chr"$chr"."$ind"_ID_length.txt > $ind/Chr"$chr"."$ind"_combined_ave_anc_pos.txt



    done
done

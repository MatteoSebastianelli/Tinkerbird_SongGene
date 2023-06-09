##########################################################

# running ELAI - example for chr. 1

# We applied 2 upper-layer clusters (-C) and 10 low-layer clusters (-c)
#Given the uncertainty related to the timing of the admixture event, 
#we investigated four possible values of the admixture generations parameter (i.e. -mg), 
#estimating ancestry scores assuming 5, 10, 15 and 20 generations since the admixture event. 
#For each migration parameter, we ran three independent runs with 30 EM steps (-s)

##########################################################

# CHROMOSOME 1
## mg5, 3 runs
elai-lin -g ref_extoni_for_elai_chr01.recode.geno.txt -p 10 -g ref_pusillus_for_elai_chr01.recode.geno.txt -p 11 -g chrysopus_to_infer_for_elai_chr01.recode.geno.txt -p 1 -pos chrysopus_to_infer_for_elai_chr01.recode.pos.txt -s 30 -o chr01_run1_mg5 -C 2 -c 10 -mg 5
elai-lin -g ref_extoni_for_elai_chr01.recode.geno.txt -p 10 -g ref_pusillus_for_elai_chr01.recode.geno.txt -p 11 -g chrysopus_to_infer_for_elai_chr01.recode.geno.txt -p 1 -pos chrysopus_to_infer_for_elai_chr01.recode.pos.txt -s 30 -o chr01_run2_mg5 -C 2 -c 10 -mg 5
#/home/vonholdt/VONHOLDT/BIN/elai/elai-lin -g ref_extoni_for_elai_chr01.recode.geno.txt -p 10 -g ref_pusillus_for_elai_chr01.recode.geno.txt -p 11 -g chrysopus_to_infer_for_elai_chr01.recode.geno.txt -p 1 -pos chrysopus_to_infer_for_elai_chr01.recode.pos.txt -s 30 -o chr01_run3_mg5 -C 2 -c 10 -mg 5

## mg10, 3 runs
elai-lin -g ref_extoni_for_elai_chr01.recode.geno.txt -p 10 -g ref_pusillus_for_elai_chr01.recode.geno.txt -p 11 -g chrysopus_to_infer_for_elai_chr01.recode.geno.txt -p 1 -pos chrysopus_to_infer_for_elai_chr01.recode.pos.txt -s 30 -o chr01_run1_mg10 -C 2 -c 10 -mg 10
elai-lin -g ref_extoni_for_elai_chr01.recode.geno.txt -p 10 -g ref_pusillus_for_elai_chr01.recode.geno.txt -p 11 -g chrysopus_to_infer_for_elai_chr01.recode.geno.txt -p 1 -pos chrysopus_to_infer_for_elai_chr01.recode.pos.txt -s 30 -o chr01_run2_mg10 -C 2 -c 10 -mg 10
elai-lin -g ref_extoni_for_elai_chr01.recode.geno.txt -p 10 -g ref_pusillus_for_elai_chr01.recode.geno.txt -p 11 -g chrysopus_to_infer_for_elai_chr01.recode.geno.txt -p 1 -pos chrysopus_to_infer_for_elai_chr01.recode.pos.txt -s 30 -o chr01_run3_mg10 -C 2 -c 10 -mg 10

## mg15, 3 runs
elai-lin -g ref_extoni_for_elai_chr01.recode.geno.txt -p 10 -g ref_pusillus_for_elai_chr01.recode.geno.txt -p 11 -g chrysopus_to_infer_for_elai_chr01.recode.geno.txt -p 1 -pos chrysopus_to_infer_for_elai_chr01.recode.pos.txt -s 30 -o chr01_run1_mg15 -C 2 -c 10 -mg 15
elai-lin -g ref_extoni_for_elai_chr01.recode.geno.txt -p 10 -g ref_pusillus_for_elai_chr01.recode.geno.txt -p 11 -g chrysopus_to_infer_for_elai_chr01.recode.geno.txt -p 1 -pos chrysopus_to_infer_for_elai_chr01.recode.pos.txt -s 30 -o chr01_run2_mg15 -C 2 -c 10 -mg 15
elai-lin -g ref_extoni_for_elai_chr01.recode.geno.txt -p 10 -g ref_pusillus_for_elai_chr01.recode.geno.txt -p 11 -g chrysopus_to_infer_for_elai_chr01.recode.geno.txt -p 1 -pos chrysopus_to_infer_for_elai_chr01.recode.pos.txt -s 30 -o chr01_run3_mg15 -C 2 -c 10 -mg 15

# mg20, 3 runs
elai-lin -g ref_extoni_for_elai_chr01.recode.geno.txt -p 10 -g ref_pusillus_for_elai_chr01.recode.geno.txt -p 11 -g chrysopus_to_infer_for_elai_chr01.recode.geno.txt -p 1 -pos chrysopus_to_infer_for_elai_chr01.recode.pos.txt -s 30 -o chr01_run1_mg20 -C 2 -c 10 -mg 20
elai-lin -g ref_extoni_for_elai_chr01.recode.geno.txt -p 10 -g ref_pusillus_for_elai_chr01.recode.geno.txt -p 11 -g chrysopus_to_infer_for_elai_chr01.recode.geno.txt -p 1 -pos chrysopus_to_infer_for_elai_chr01.recode.pos.txt -s 30 -o chr01_run2_mg20 -C 2 -c 10 -mg 20
elai-lin -g ref_extoni_for_elai_chr01.recode.geno.txt -p 10 -g ref_pusillus_for_elai_chr01.recode.geno.txt -p 11 -g chrysopus_to_infer_for_elai_chr01.recode.geno.txt -p 1 -pos chrysopus_to_infer_for_elai_chr01.recode.pos.txt -s 30 -o chr01_run3_mg20 -C 2 -c 10 -mg 20

# The main output is stored in .ps21.txt files, which I then move in a folder specific for
#each migration parameter and run (e.g. mg05_run1, mg05_run2, etc...)

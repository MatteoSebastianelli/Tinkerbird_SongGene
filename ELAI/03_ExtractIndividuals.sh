##########################################################

# Extracting individual information from ELAI output - example for first 3 individuals

# In the .ps21.txt files, there are as many rows as the number of individuals, 
#these are ordered as in the initial PLINK files. We want to extract information for each
#individual per run per migration parameter.
#We then move all the files in the same folder (i.e. runs)

##########################################################


cd /mg05_run1 #For this task, I run 12 jobs, one for every directory containing ELAI output
#cd /mg05_run2
#cd /mg05_run3

#cd /mg10_run1
#cd /mg10_run2
#cd /mg10_run3

#cd /mg15_run1
#cd /mg15_run2
#cd /mg15_run3

#cd /mg20_run1
#cd /mg20_run2
#cd /mg20_run3


#AR93112
files=*ps21.txt
for f in $files
do
pre=$(echo $f |awk -F'.' '{print $1 }')
suf=".AR93112.txt" # Add extension to file name
sed -n 1p $f > $pre$suf # Extract 1st row
done

mv *.AR93112.txt /move/to/common/folder/runs

#AR93113
files=*ps21.txt
for f in $files
do
pre=$(echo $f |awk -F'.' '{print $1 }')
suf=".AR93113.txt"
sed -n 2p $f > $pre$suf # Extract 2nd row
done

mv *.AR93113.txt /move/to/common/folder/runs


#AR93114
files=*ps21.txt
for f in $files
do
pre=$(echo $f |awk -F'.' '{print $1 }')
suf=".AR93114.txt"
sed -n 3p $f > $pre$suf # Extract 3rd row
done

mv *.AR93114.txt /move/to/common/folder/runs

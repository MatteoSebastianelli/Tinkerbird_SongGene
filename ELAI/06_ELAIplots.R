require(ggplot2)
require(tidyverse)

##### Load reference genome ####
tink<- read.table("/tigress/VONHOLDT/matteo/ref_genome/bPogPus1.pri.cur.20200514.fastadoc.txt",
                  header = F) 

colnames(tink)<-c("Chromosome", "Length") #name the first two columns
head(tink)

tink1<-tink[,c(1,2)] # drop columns 3 to 5
tink1$Chromosome

#### Extract chromosomes only - no scaffolds ####
tink_genome<-subset(tink1, tink1$Chromosome=="SUPER_1" | 
                      tink1$Chromosome=="SUPER_2" |
                      tink1$Chromosome=="SUPER_3" |
                      tink1$Chromosome=="SUPER_4" |
                      tink1$Chromosome=="SUPER_5" |
                      tink1$Chromosome=="SUPER_6" |
                      tink1$Chromosome=="SUPER_7" |
                      tink1$Chromosome=="SUPER_8" |
                      tink1$Chromosome=="SUPER_9" |
                      tink1$Chromosome=="SUPER_10" |
                      tink1$Chromosome=="SUPER_11" |
                      tink1$Chromosome=="SUPER_12" |
                      tink1$Chromosome=="SUPER_13" |
                      tink1$Chromosome=="SUPER_14" |
                      tink1$Chromosome=="SUPER_15" |
                      tink1$Chromosome=="SUPER_16" |
                      tink1$Chromosome=="SUPER_17" |
                      tink1$Chromosome=="SUPER_18" |
                      tink1$Chromosome=="SUPER_19" |
                      tink1$Chromosome=="SUPER_20" |
                      tink1$Chromosome=="SUPER_21" |
                      tink1$Chromosome=="SUPER_22" |
                      tink1$Chromosome=="SUPER_23" |
                      tink1$Chromosome=="SUPER_24" |
                      tink1$Chromosome=="SUPER_25" |
                      tink1$Chromosome=="SUPER_26" |
                      tink1$Chromosome=="SUPER_27" |
                      tink1$Chromosome=="SUPER_28" |
                      tink1$Chromosome=="SUPER_29" |
                      tink1$Chromosome=="SUPER_30" |
                      tink1$Chromosome=="SUPER_31" |
                      tink1$Chromosome=="SUPER_32" |
                      tink1$Chromosome=="SUPER_33" |
                      tink1$Chromosome=="SUPER_34" |
                      tink1$Chromosome=="SUPER_35" |
                      tink1$Chromosome=="SUPER_36" |
                      tink1$Chromosome=="SUPER_37" |
                      tink1$Chromosome=="SUPER_38" |
                      tink1$Chromosome=="SUPER_39" |
                      tink1$Chromosome=="SUPER_40" |
                      tink1$Chromosome=="SUPER_41" |
                      tink1$Chromosome=="SUPER_42" |
                      tink1$Chromosome=="SUPER_43" |
                      tink1$Chromosome=="SUPER_44" |
                      tink1$Chromosome=="SUPER_W" |
                      tink1$Chromosome=="SUPER_Z")

# Rename chromosomes
tink_genome$Chromosome<-factor(tink_genome$Chromosome, levels = c("SUPER_1","SUPER_2", "SUPER_3", "SUPER_4", "SUPER_5", "SUPER_6", "SUPER_7", "SUPER_8", "SUPER_9", "SUPER_10",
                                                                  "SUPER_11", "SUPER_12", "SUPER_13", "SUPER_14", "SUPER_15", "SUPER_16", "SUPER_17", "SUPER_18", "SUPER_19", "SUPER_20",
                                                                  "SUPER_21", "SUPER_22", "SUPER_23", "SUPER_24", "SUPER_25", "SUPER_26", "SUPER_27", "SUPER_28", "SUPER_29", "SUPER_30",
                                                                  "SUPER_31", "SUPER_32", "SUPER_33", "SUPER_34", "SUPER_35", "SUPER_36", "SUPER_37", "SUPER_38", "SUPER_39", "SUPER_40",
                                                                  "SUPER_41", "SUPER_42", "SUPER_43", "SUPER_44","SUPER_W", "SUPER_Z"), 
                               labels = c("1","2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                                          "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39",
                                          "40", "41", "42", "43", "44", "45", "46"))

tink_genome <- tink_genome[order(tink_genome$Chromosome),] #order by chromosome
str(tink_genome)

tink_genome1<-subset(tink_genome, tink_genome$Chromosome!="45") # Remove chr W (i.e.45)

par(lwd = 1)
tink_genome$Length_mb<-(tink_genome$Length)/1000000
tink_genome1$Length_mb<-(tink_genome1$Length)/1000000

genome_plot <-barplot(tink_genome1$Length_mb, col="grey80", border=1, lwd=2,
                      ylab="", xlab="", las=1)

#### Upoload individual genomes ####
d <- read.table("""$ind"_AllChr_combined_ave_anc_pos.txt", # This file was produced in the previous script (i.e. 05)
                fill = T,
                header = F)

colnames(d) <- c("CHR", "POS", "ANC", "ID")

d$POS1<-(d$POS)/1000000

head(d)

d1<-subset(d, d$CHR!="45")

d1$CHR<-as.factor(d1$CHR)

col=ifelse(d1$ANC >= 0.5 & d1$ANC < 1.5, 'gray80',
           ifelse (d1$ANC >=1.5, 'darkred', '#FEB24C'))
##CC0000 -> red
#FEB24C -> yellow
#FFFFCC -> beige

d1$CHR<-as.numeric(d1$CHR) #convert chromosome as numeric

#add over original plot and save
pdf("ELAI_plot.pdf", paper = "a4r")

genome_plot <-barplot(tink_genome1$Length_mb, col="grey80", border=1, lwd=2,
                      ylab="", xlab="", las=1)

with(d1,
     segments(
       genome_plot[CHR,]-0.4,
       POS1,
       genome_plot[CHR,]+0.4,
       POS1,
       col=col,
       lwd=1, 
       lend=1
     )
)

# Close the pdf file
dev.off() 


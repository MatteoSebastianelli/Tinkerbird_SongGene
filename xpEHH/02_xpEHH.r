##########################################################

# Compute 
# This script follows Mark Ravinet's speciation genomics tutorial (https://speciationgenomics.github.io/haplotypes/)

##########################################################

#R version 4.1.3 (2022-03-10) -- "One Push-Up"

require(rehh)
require(tidyr)
require(dplyr)
require(ggplot2)
require(patchwork)
require(RAINBOWR)
require(qvalue)
require(writexl)
require(readxl)

#### Import data ####
# Import allopatric YFT and RFT separately
yft <- data2haplohh(hap_file = "../xpEHH/Chr25_yft.vcf.gz",
                    polarize_vcf = F)

rft <- data2haplohh(hap_file = "../xpEHH/Chr25_rft.vcf.gz",
                    polarize_vcf = F)
                    
#### Performing a haplotype genome scan ####
yft_scan <- scan_hh(yft, polarized = FALSE)
rft_scan <- scan_hh(rft, polarized = FALSE)

#### Calculating xpEHH ####
chrysopus_xpehh <- ies2xpehh(yft_scan, rft_scan,
                       popname1 = "chrysoconus", popname2 = "pusillus",
                       include_freq = T)
                       
#### Highlight significant snps from gemma ####
snpsOfInterest <- c(5022427, 6851380, 6854983, 6957543, 7737388, 8886686, 
                                  9449012, 9589214, 9602675, 9664989, 9791236, 9797804,
                                  10249778, 11568784, 17821223)

names(snpsOfInterest) <- c('rs')
str(snpsOfInterest)

chrysopus_xpehh <- chrysopus_xpehh %>% 
  mutate(is_highlight=ifelse(POSITION %in% snpsOfInterest, "yes", "no"))

any(chrysopus_xpehh$is_highlight == "yes")

str(chrysopus_xpehh)

a <- subset(chrysopus_xpehh, chrysopus_xpehh$is_highlight=="yes")

# get 2.5% and 97.5% quantiles
quantile(chrysopus_xpehh$XPEHH_chrysoconus_pusillus, probs = c(0.025, 0.5, 0.975), na.rm = T)

#### Plot ####
chrysopus_xpehh$pos_Mb <- chrysopus_xpehh$POSITION/1000000

ggplot(chrysopus_xpehh, aes(pos_Mb, XPEHH_chrysoconus_pusillus)) + 
  geom_point(col = "gray75", alpha = 0.7, shape = 1) +
  geom_point(data=subset(chrysopus_xpehh, is_highlight=="yes"), color="darkred", size=2, shape = 19) +
  geom_hline(yintercept = -2.01129124, linetype = "dashed", size = 1.2, col = "black") +
  geom_hline(yintercept = 1.85729582, linetype = "dashed", size = 1.2, col = "black") +
  theme_classic() +
  labs(x = "Chromosome 25", y = "xp-EHH") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(colour = "black"),
        legend.position = "none")
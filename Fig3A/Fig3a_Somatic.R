###
## Association between MSI and germline/somatic mutations
## 
###
library(data.table)
library(tidyverse)
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MSI_MMR/MMR_MSI/")

#all = fread("All_Pathogenic_variants.tsv") #Read variants
d = filter(all, state == "Somatic")

#Filter variants by VAF and gnomad populaton frequency
d = filter(d, VAF>0.20)
d = filter(d, Gnomad_AF < 0.005|is.na(Gnomad_AF))

#Sum over cancertypes and calculate the proportion of MSI in each gene in each study
d2 = d%>%
  group_by(Study, Gene, msiStatus)%>%
  mutate(
    n = n_distinct(Sample_ID)
  )%>%
  distinct(Study, Gene, msiStatus, .keep_all = T)%>%
  pivot_wider(names_from = msiStatus, values_from = n)%>%
  group_by(Study, Gene)%>%
  mutate(
    n_MSI = sum(MSI, na.rm = T),
    n_MSS = sum(MSS, na.rm = T),
    n_mutated = sum(n_MSI, n_MSS)/2,
    prop_MSI = n_MSI/n_mutated,
    max_hit = max(n_MSI, n_MSS)
  )%>%
  ungroup()%>%
  distinct(Study, Gene, .keep_all = T)%>%
  dplyr::select(-cancertype, -MSI, -MSS)%>%
  group_by(Study)%>%
  filter(n_mutated >= 3)%>%
  distinct()%>%
  mutate(
    median_prop = median(prop_MSI)
  )

table(d2$median_prop) #See the medians

#Select genes for foreground..
goi = readxl::read_excel("Data/Supplementary tables.xlsx", sheet = 7)
goi = goi[goi$`Core gene`== "Yes",]

d3 = filter(d2, Gene %in% goi$Gene)%>%
  filter(n_MSI>0, n_mutated >= 1)


#Binom-test per sample
for(i in 1:nrow(d3)){
  ti = binom.test(x = d3$n_MSI[i], n = d3$n_mutated[i], p = d3$median_prop[i], alternative = "greater")
  d3$"p_val"[i] = ti$p.value
}

#Combine p-values using fishers method (sumlog). 
source("Fig3A/sumlog.R")
d4 = dplyr::select(d3, Gene, Study, n_MSI, n_mutated, prop_MSI, p_val)%>%
  group_by(Gene)%>%
  mutate(
    p_val2 = p_val,
    p_val = ifelse(n_distinct(p_val)>1,
                   sumlog(c(p_val))$p, 
                   sum(p_val, na.rm = T)),
    p_val = ifelse(p_val == 2, 1, p_val)
  )%>%
  pivot_wider(names_from = Study, values_from = c(n_MSI, n_mutated, prop_MSI, p_val2))

#Get a shared count across the two data frames
d4 = d4%>%
  mutate(
    sum_n_MSI = n_MSI_PCAWG + n_MSI_HMF,
    sum_n_mutated = n_mutated_PCAWG + n_mutated_HMF,
    sum_prop = sum_n_MSI/sum_n_mutated
  )

#FDR correction
d4$q_val = p.adjust(d4$p_val, method = "fdr", n = nrow(d4)) 

#Write table to combine with germline version for plotting
write.table(d4, "Figures/Fig3_new/somatic_fig3.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

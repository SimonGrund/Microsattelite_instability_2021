library(data.table)
library(tidyverse)
library(scales)
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MSI_MMR/Github/Microsattelite_instability_2021")

#Load MSI status from the supplementary table, sheet 3
d = readxl::read_excel("Data/Supplementary tables.xlsx", sheet = 3)
colnames(d) = c("Donor_ID", "Sample_ID", "cancertype", "msiStatus", "S.ind")

#Calculate the proportion of MSI in each cancertype and use a binomial test to compare to the global proportion (183/6098)
d2 = d%>%
  group_by(cancertype)%>%
  mutate(
    n_ct = n_distinct(Sample_ID),
    n = sum(msiStatus == "MSI")
  )%>%
  ungroup()%>%
  distinct(cancertype, n_ct, n)%>%
  mutate(prop = n/n_ct)%>%
  distinct(cancertype, .keep_all = T)%>%
  rowwise()%>%
  mutate(
    p_val = binom.test(x = n, n = n_ct, p = 183/6098, alternative = "greater")$p.value,
  )

d3 = filter(d2, n >0)
d3$q_val = p.adjust(d3$p_val, method = "fdr")

d3 = mutate(d3, 
            sign = ifelse(q_val < 0.10, "q \U2264 0.10", "NS"),
            prop = n/n_ct
)

gA = ggplot(d3, aes(x = reorder(paste(cancertype," (MSI: ", n, " / ", n_ct, ")", sep =""), -prop), y = prop))+
  geom_histogram(stat = "identity",aes(fill = sign), color = "black", width = 0.6, show.legend = T, size = 0.2)+
  scale_fill_manual(values = c("darkblue", "grey60"),
                    breaks = c("q \U2264 0.10", "NS"))+
  geom_text(data = filter(d3, q_val < 0.01), aes(label = scales::scientific(q_val, digits = 2,prefix = " q \U2264 ")),
            vjust = 0.5, hjust = 0, angle = 90, nudge_y = 0, size = 2.5, family = "Arial")+
  geom_text(data = filter(d3, q_val > 0.01), aes(label = paste(" q \U2264 ", round(q_val, 2), sep = "")),
            vjust = 0.5, hjust = 0, angle = 90, nudge_y = 0, size = 2.5, family = "Arial")+
  theme_bw(base_size = 9, base_family = "Arial") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.position = c(0.8, 0.8),
        legend.direction = "horizontal",
        plot.margin = margin(5,5,5,50),
        axis.line = element_line(color = NA))+
  xlab("")+
  ylab("Proportion of samples with MSI")+
  scale_y_continuous(limits =c(0,0.42))+
  guides(fill = guide_legend(title = ""))
gA

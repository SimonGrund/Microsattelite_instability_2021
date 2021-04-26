library(data.table)
library(tidyverse)
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MSI_MMR/MMR_MSI/")

ds = fread("Figures/Fig3_new/somatic_fig3.tsv")
dg = fread("Figures/Fig3_new/germline_fig3.tsv")

ds$"state" = "Somatic"
dg$"state" = "Germline"
d = bind_rows(ds,dg)

d = mutate(d, sign = ifelse(q_val < 0.1, "q \U2264 0.10", "NS"))

ggplot(d, aes(x = sum_prop, y = -log(p_val), 
              label = paste(Gene,": ", sum_n_MSI, " / ", sum_n_mutated, sep = ""),
       col = sign))+
  geom_point(size = 0.5, alpha = 0.5)+
 ggrepel::geom_text_repel(data = filter(d, q_val<0.1), size = 2.5, family = "Arial", segment.size = 0.1,
                          nudge_y = 1, min.segment.length = 0.001, col = "darkblue")+
  # geom_text(data = filter(d, q_val>0.1),
  #                          aes(label = Gene),
  #           #segment.size = 0.1,
  #                          size = 2.5, family = "Arial", col = "grey20",
  #          nudge_y = -0.1)+
  ggrepel::geom_text_repel(data = filter(d, q_val>0.1), aes(label = Gene), size = 2, family = "Arial", segment.size = 0.1,
                           nudge_y = 0.1, min.segment.length = 0.01)+
  theme_bw(base_size = 7, base_family = "Arial")+
  theme(
    axis.text.x = element_text(size = 7),
    axis.title = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    plot.margin = margin(5,5,5,50),
    axis.line = element_line(color = NA),
   strip.text = element_text(size = 7),
   panel.grid.major = element_blank(),
   panel.grid.minor = element_blank()
  )+
  xlab("Proportion of mutated samples with MSI")+
  ylab("-log(p-value)")+
  facet_wrap(~state, nrow = 2, scales = "free_x")+
#  facet_grid(rows = vars(state), scales = "free_x")+
  scale_color_manual(values = c("grey20", "darkblue"),
                    breaks = c("NS", "q \U2264 0.10"))

ggsave("Figures/Fig3_new/Fig3a.pdf", device = cairo_pdf, width = 8, height = 5)

library(flowCore)
library(flowWorkspace)
library(ggcyto)
library(flowAI)
library(cytotidyr) # this is supposed to help.
library(CytoML)
library(dplyr)
library(tidyr)
library(ggplot2)
library(uwot)
library(Rtsne)
library(scales)
library(flowViz)
library(cytotidyr) # this allows us to convert a flowset into a dataframe
library(ggthemes)
library(gridExtra)

# C-c RET is now bound as the ess-eval-region-or-line-and-step function.

cs <- c()
try(cs <- read.flowSet(path = "./resultsQC/", pattern = ".fcs"))
if (length(cs) == 0) {
  theme_set(theme_clean())
  read.flowSet(path = "../20200810NJ-ACEFlow/good/", pattern = ".fcs") %>%
    flow_auto_qc() %>% 
    flowSet_to_cytoset() -> cs
}
new_sample_names <- c()
for (i in 1:length(cs)) {
  new_sample_names <- c(new_sample_names, 
                        keyword(cs[[i]])$`TUBE NAME`)
}
sampleNames(cs) <- new_sample_names
gs <- GatingSet(cs)

pData(gs) <- cbind(pData(gs),row.names(pData(gs)))
colnames(pData(gs))[2] <- "tube names"

# # Figure out your gate locations
# fs <- gs_pop_get_data(gs, "scatter")
# c_x <- "SSC-W"
# c_y <- "SSC-H"
# p_x <- autoplot(fs, x = c_x)
# p_y <- autoplot(fs, x = c_y)
# p_xy <- autoplot(fs, x = c_x, y = c_y)
# g_x <- 1.7e+05
# g_y <- 6e+04
# p_x + geom_vline(xintercept = g_x)
# p_y + geom_vline(xintercept = g_y)
# p_xy + geom_vline(xintercept = g_x) + geom_hline(yintercept = g_y) + ggcyto_par_set(limits = "instrument")

scatter <- matrix(byrow = TRUE,
               ncol = 2,
               data = c(
                 5e+04, 1.5e+04,
                 5e+04, 5.0e+04,
                 1e+05, 1.2e+05,
                 2e+05, 1.2e+05,
                 2e+05, 1.5e+04
               ))
colnames(scatter) <- c("FSC-A","SSC-A")
scatter <- polygonGate(scatter, filterId = "scatter")
gs_pop_add(gs, scatter, name = "scatter", parent = "root")

ggcyto(gs, 
       aes(`FSC-A`, `SSC-A`), 
       subset = "root", 
       filter = marginalFilter) + 
  geom_bin2d(bins = 256) + 
  geom_gate("scatter") + 
  labs(title = "Root Scatter Plots.pdf") +
  facet_wrap(~`tube names`)
  

fsc_singlet <- rectangleGate(filterId = "FSC-singlet", 
                             "FSC-W" = c(9e+04, 1.7e+05),
                             "FSC-H" = c(3e+04, 7e+04)
                             )
gs_pop_add(gs, 
           fsc_singlet, 
           name = "FSC-singlet", 
           parent = "scatter")

ssc_singlet <- rectangleGate(filterId = "SSC-singlet", 
                             "SSC-W" = c(0.9e+05, 1.7e+05),
                             "SSC-H" = c(10, 6e+04)
)
gs_pop_add(gs, 
           ssc_singlet, 
           name = "SSC-singlet", 
           parent = "FSC-singlet")
recompute(gs)

gs_pop_get_data(gs[["US"]], "SSC-singlet") %>% 
  fortify() %>% 
  pull(`FITC-A`) %>%
  quantile(.99) -> max_FITC

mace_bound <- rectangleGate(filterId="MACE-bound","FITC-A"=c(max_FITC,Inf))
gs_pop_add(gs, 
           mace_bound, 
           name = "MACE-bound", 
           parent = "SSC-singlet")
recompute(gs)

nodes <- gs_get_pop_paths(gs, path = "auto")[2:4]

pdf("Root Scatter Plots.pdf", 7, 7)
ggcyto(gs, 
       aes(`FSC-A`, `SSC-A`), 
       subset = "root", 
       filter = marginalFilter) + 
  geom_bin2d(bins = 256) + 
  ggcyto_par_set(limits = "instrument") +
  geom_gate("scatter") + 
  geom_stats("scatter", type = "percent", adjust = 0.8) +
  labs(title = "Root Scatter Plots") +
  facet_wrap(~`tube names`)
dev.off()
pdf("FSC-singlet Plots.pdf", 7, 7)
ggcyto(gs, 
       aes(`FSC-W`, `FSC-H`), 
       subset = "scatter", 
       filter = marginalFilter) + 
  geom_bin2d(bins = 256) + 
  ggcyto_par_set(limits = "instrument") +
  geom_gate("FSC-singlet") + 
  geom_stats("FSC-singlet", type = "percent", adjust = 0.1) +
  labs(title = "FSC-singlet Plots") +
  facet_wrap(~`tube names`)
dev.off()
pdf("SSC-singlet Plots.pdf", 7, 7)
ggcyto(gs, 
       aes(`SSC-W`, `SSC-H`), 
       subset = "FSC-singlet", 
       filter = marginalFilter) + 
  geom_bin2d(bins = 256) + 
  ggcyto_par_set(limits = "instrument") +
  geom_gate("SSC-singlet") + 
  geom_stats("SSC-singlet", type = "percent", adjust = 0.1) +
  labs(title = "SSC-singlet Plots") +
  facet_wrap(~`tube names`)
dev.off()
pdf("MACE-binding Plots.pdf", 7, 7)
ggcyto(gs, 
       aes(`FITC-A`, `SSC-A`), 
       subset = "SSC-singlet") + 
  geom_bin2d(bins = 256) + 
  scale_x_flowjo_biexp(neg = 1, widthBasis = -20) +
  geom_gate("MACE-bound") + 
  geom_stats("MACE-bound", type = "percent", adjust = 0.1) +
  labs(title = "MACE-binding Plots", subtitle = "Gated at the 99th percentile of the unstained sample") +
  facet_wrap(~`tube names`)
dev.off()

ggcyto(gs, 
       aes(`FITC-A`, `SSC-A`), 
       subset = "SSC-singlet") + 
  geom_bin2d(bins = 256) + 
  scale_x_flowjo_biexp(neg = 1, widthBasis = -20) +
  geom_gate("MACE-bound") + 
  geom_stats("MACE-bound", type = "percent", adjust = 0.1) +
  labs(title = "MACE-binding", subtitle = "Gated at the 99th percentile of the unstained sample") +
  facet_wrap(~`tube names`)
ggsave("MACE-binding.pdf", width = 5, height = 5)

# Root Scatter Plots
# FSC-singlet Plots
# SSC-singlet Plots
# Mace-binding Plots

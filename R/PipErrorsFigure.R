library(dplyr)
library(ggplot2)
library(limma)
library(tibble)
library(matrixStats)
library(scales)
library(patchwork)
library(cowplot)
library(ggpubr)
library(reshape2)
library(stringr)
library(plyr)

# Read in the output from FDPAnalysis.ipynb ----
results_path <- r"(...\FDP_Analysis_Results.tsv)"
results <- read.csv(results_path, sep = '\t')

#subset the relevant portion of the output
results_mini <- results[c("software", "fper", "nper", "fder")]
results_melt <- melt(results_mini, value.name = "error_rate", columns = c("fper", "nper", "fder"))


#Work up the data (Figure 3A)----
results_melt$Dataset <- strsplit(results_melt$software, "_") %>% lapply(function(x) x[1]) %>% unlist()
results_melt <- results_melt[order(results_melt$Dataset),]
results_melt$Dataset[results_melt$Dataset == "Gygi"] <- "Lim et al."
results_melt$Dataset[results_melt$Dataset == "Inhouse"] <- "E. coli"
results_melt$Dataset[results_melt$Dataset == "Kelly"] <- "Single Cell"

results_melt$Dataset <- factor(results_melt$Dataset, levels = c("Lim et al.", "Single Cell", "E. coli"))

results_melt$variable <- revalue(results_melt$variable, c("fper" = "Foreign Peak Errors", "nper" = "Native Peak Errors", 
                                          "fder" = "False Detection Errors"))
results_melt$variable <- factor(results_melt$variable, levels = rev(c("Foreign Peak Errors", 
                                                      "Native Peak Errors",
                                                      "False Detection Errors")))


results_melt$SearchEngine <- strsplit(results_melt$software, "_") %>% lapply(function(x) x[-1]) %>%
  lapply(function(x) paste(x, collapse="_")) %>% unlist()

results_melt$PipFdr <-  strsplit(results_melt$SearchEngine, "_") %>% lapply(function(x) tail(x, n=1)) %>% unlist()
results_melt$PipFdr[results_melt$PipFdr == "v1"] <- "100"
results_melt$PipFdr[results_melt$PipFdr == "MaxQuant"] <- "100"
results_melt$PipFdr[results_melt$PipFdr == "100"] <- "N/A"
results_melt$PipFdr <- results_melt$PipFdr %>% factor()

results_melt$SearchEngine[grepl("Flash_v2", results_melt$SearchEngine)] <- "FlashLFQ +\nPIP-ECHO"
results_melt$SearchEngine[grepl("Flash_v1", results_melt$SearchEngine)] <- "Flash\nv1.0"
results_melt$SearchEngine[grepl("IonQuant", results_melt$SearchEngine)] <- "IonQuant"
results_melt$SearchEngine[grepl("MaxQuant", results_melt$SearchEngine)] <- "Max-\nQuant"

results_melt$SearchEngine <- factor(results_melt$SearchEngine, 
                            levels = c("FlashLFQ +\nPIP-ECHO",  "IonQuant", "Max-\nQuant", "Flash\nv1.0"))


# Plot and combine the stacked bar plots for each dataset ----
gygi <- results_melt[results_melt$Dataset == "Lim et al.", ]
sc <- results_melt[results_melt$Dataset == "Single Cell", ]
inhouse <- results_melt[results_melt$Dataset == "E. coli", ]

gygi_fdp_leg <- ggplot(gygi, aes(fill = variable, y = error_rate, x = PipFdr)) +
  geom_bar(position="stack", stat = "identity") +
  facet_grid( ~ SearchEngine, space = "free", scales = "free") +
  coord_cartesian(ylim = c(0, 6.5)) +
  ylab("Percent False Discoveries") +
  xlab("PIP FDR") + 
  theme_bw() +
  theme(panel.grid.major.y = element_line(color = "darkgray")) +
  scale_y_continuous(breaks = c(1, 2.5, 5), labels = c(1, 2.5, 5)) +
  scale_fill_manual(values = alpha(rev(c('#D55E00', '#E69F00','#F0E442')), 0.8), 
                    name = "Error Type")

fdp_legend <- get_legend(gygi_fdp_leg)


ax_text_size = 15
ax_title_size = 20
plot_title_size = 24
strip_text_size = 13

y_start = 0
y_end = 6

gygi_fdp <- ggplot(gygi, aes(fill = variable, y = error_rate, x = PipFdr)) +
  geom_bar(position="stack", stat = "identity") +
  facet_grid( ~ SearchEngine, space = "free", scales = "free") +
  coord_cartesian(ylim = c(y_start, y_end)) +
  ylab("") +
  xlab("PIP FDR") + 
  theme_bw() +
  theme(panel.grid.major.y = element_line(color = "darkgray"),
        legend.position = "none",
        axis.text = element_text(size=ax_text_size),
        axis.title=element_text(size=ax_title_size),
        plot.title = element_text(size=plot_title_size),
        strip.text.x = element_text(size = strip_text_size)) +
  scale_y_continuous(breaks = c(1, 2.5, 5), labels = c(1, 2.5, 5)) +
  scale_fill_manual(values = alpha(rev(c('#D55E00', '#E69F00','#F0E442')), 0.8), name = "Error Type") +
  ggtitle("Lim et al.")

sc_fdp <- ggplot(sc, aes(fill = variable, y = error_rate, x = PipFdr)) +
  geom_bar(position="stack", stat = "identity") +
  facet_grid( ~ SearchEngine, space = "free", scales = "free") +
  coord_cartesian(ylim = c(y_start, y_end)) +
  ylab("Percent False Discoveries") +
  xlab("PIP FDR") + 
  theme_bw() +
  theme(panel.grid.major.y = element_line(color = "darkgray"),
        legend.position = "none",
        axis.text = element_text(size=ax_text_size),
        axis.title=element_text(size=ax_title_size),
        plot.title = element_text(size=plot_title_size),
        strip.text.x = element_text(size = strip_text_size)) +
  scale_y_continuous(breaks = c(1, 2.5, 5), labels = c(1, 2.5, 5)) +
  scale_fill_manual(values = alpha(rev(c('#D55E00', '#E69F00','#F0E442')), 0.8), name = "Error Type") +
  ggtitle("Single cell")

inhouse_fdp <- ggplot(inhouse, aes(fill = variable, y = error_rate, x = PipFdr)) +
  geom_bar(position="stack", stat = "identity") +
  facet_grid( ~ SearchEngine, space = "free", scales = "free") +
  coord_cartesian(ylim = c(y_start, y_end)) +
  ylab("") +
  xlab("PIP FDR") + 
  theme_bw() +
  theme(panel.grid.major.y = element_line(color = "darkgray"),
        legend.position = "none",
        axis.text = element_text(size=ax_text_size),
        axis.title=element_text(size=ax_title_size),
        plot.title = element_text(size=plot_title_size),
        strip.text.x = element_text(size = strip_text_size)) +
  scale_y_continuous(breaks = c(1, 2.5, 5), labels = c(1, 2.5, 5)) +
  scale_fill_manual(values = alpha(rev(c('#D55E00', '#E69F00','#F0E442')), 0.8), name = "Error Type") +
  ggtitle("E. coli")

fdp_grid <- plot_grid( sc_fdp, inhouse_fdp, gygi_fdp, ncol = 3)
fdp_grid

out_dir <- r"(C:\...\)"

ggsave(filename = paste0(out_dir, "fdp_plot_7772", ".png"),
       plot = fdp_grid,
       units = "px",
       width = 5000,
       height = 2000)

# Count the number of PIP Events for each software (Figure 3C) ----

pip_counts <- results[c("software", "total")]

pip_counts$Dataset <- strsplit(pip_counts$software, "_") %>% lapply(function(x) x[1]) %>% unlist()
pip_counts <- pip_counts[order(pip_counts$Dataset),]
pip_counts$Dataset[pip_counts$Dataset == "Gygi"] <- "Lim et al."
pip_counts$Dataset[pip_counts$Dataset == "Inhouse"] <- "E. coli"
pip_counts$Dataset[pip_counts$Dataset == "Kelly"] <- "Single Cell"

pip_counts$Dataset <- factor(pip_counts$Dataset, levels = c("Lim et al.", "Single Cell", "E. coli"))


pip_counts$SearchEngine <- strsplit(pip_counts$software, "_") %>% lapply(function(x) x[-1]) %>%
  lapply(function(x) paste(x, collapse="_")) %>% unlist()

pip_counts$PipFdr <-  strsplit(pip_counts$SearchEngine, "_") %>% lapply(function(x) tail(x, n=1)) %>% unlist()
pip_counts$PipFdr[pip_counts$PipFdr == "v1"] <- "100"
pip_counts$PipFdr[pip_counts$PipFdr == "MaxQuant"] <- "100"
pip_counts$PipFdr[pip_counts$PipFdr == "100"] <- "N/A"
pip_counts$PipFdr <- pip_counts$PipFdr %>% factor()

pip_counts$SearchEngine[grepl("Flash_v2", pip_counts$SearchEngine)] <- "FlashLFQ +\nPIP-ECHO"
pip_counts$SearchEngine[grepl("Flash_v1", pip_counts$SearchEngine)] <- "Flash\nv1.0"
pip_counts$SearchEngine[grepl("IonQuant", pip_counts$SearchEngine)] <- "IonQuant"
pip_counts$SearchEngine[grepl("MaxQuant", pip_counts$SearchEngine)] <- "Max-\nQuant"

pip_counts$SearchEngine <- factor(pip_counts$SearchEngine, 
                                  levels = c("FlashLFQ +\nPIP-ECHO",  "IonQuant", "Max-\nQuant", "Flash\nv1.0"))

gygi_pip <- pip_counts[pip_counts$Dataset == "Lim et al.", ]
sc_pip <- pip_counts[pip_counts$Dataset == "Single Cell", ]
inhouse_pip <- pip_counts[pip_counts$Dataset == "E. coli", ]



plot_pip_counts <- function(data, title = "x", y_title = "")
{
  ggplot(data, aes(y = total, x = PipFdr)) +
    geom_bar(stat = "identity", fill = alpha("#0072B2", .75)) +
    facet_grid( ~ SearchEngine, space = "free", scales = "free") +
    ylab(y_title) +
    xlab("PIP FDR") + 
    theme_bw() +
    ggtitle(title) +
    scale_y_continuous(label = comma) +
    #coord_cartesian(ylim = c(0, 300000)) +
    theme(panel.grid.major.y = element_line(color = "darkgray"),
          legend.position = "none",
          axis.text = element_text(size=ax_text_size),
          axis.title=element_text(size=ax_title_size),
          plot.title = element_text(size=plot_title_size),
          strip.text.x = element_text(size = strip_text_size))
}


gygi_pip_plot <- plot_pip_counts(gygi_pip, "Lim et al.")

sc_pip_plot <- plot_pip_counts(sc_pip, "Single cell", y_title = "Total PIP Events")

inhouse_pip_plot <- plot_pip_counts(inhouse_pip, "E. coli")

count_grid <- plot_grid(sc_pip_plot, inhouse_pip_plot, gygi_pip_plot, ncol = 3)
count_grid


out_dir <- r"(...)"

ggsave(filename = paste0(out_dir, "pip_count_plot_7772", ".png"),
       plot = count_grid,
       units = "px",
       width = 5000,
       height = 2000)

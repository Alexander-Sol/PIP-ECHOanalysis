library(dplyr)
library(ggplot2)
library(limma)
library(tibble)
library(matrixStats)
library(scales)
library(patchwork)
library(cowplot)
library(ggpubr)

# Specific Comparison IonQuant ----

# Paths
iq_pip <- "D:/PXD003881_IonStar_SpikeIn/IonQuant_Normed_Pip5/combined_ion.tsv"
iq_no_pip <- "D:/PXD003881_IonStar_SpikeIn/IonQuant_Normed_NoPip/combined_ion.tsv"

# Special functions
count_peps_at_fdp_iq <- function(peptide_table, conditions = c("A_", "E_"), fc = 3, fdp = 0.05, order_by_fc = F)
{
  cols <- colnames(peptide_table)[grepl(paste(conditions, collapse = "|"), colnames(peptide_table))]
  intensity_cols = cols[grepl("Intensity",cols)]
  i_cols_a <- intensity_cols[grepl(conditions[1], intensity_cols)]
  i_cols_e <- intensity_cols[grepl(conditions[2], intensity_cols)]

  peptide_table[peptide_table == 0] <- NA
  peptide_table["A_median"] <- apply(peptide_table[i_cols_a], MARGIN = 1, function(x) mean(log2(x), na.rm=TRUE))
  peptide_table["E_median"] <- apply(peptide_table[i_cols_e], MARGIN = 1, function(x) mean(log2(x), na.rm=TRUE))
  peptide_table["Log2_fc"] <- peptide_table["E_median"] - peptide_table["A_median"]
  peptide_table <- peptide_table[!is.na(peptide_table$Log2_fc),]
  peptide_table["Log2_fc_delta"] <- abs(peptide_table$Log2_fc - log2(fc))
  
  # remove duplicates (smae peptide, different charges)
  peptide_table$quant_count <- apply(peptide_table, 1, function(x) length(which(x %in% c("MBR", "MS/MS"))))
  peptide_table <- peptide_table %>% 
    arrange(desc(quant_count)) %>%
    group_by(Modified.Sequence) %>%
    slice_head(n=1)
  
  if(order_by_fc)
  {
    peptide_table <- peptide_table[order(peptide_table$Log2_fc, decreasing = TRUE),]
  } else
  {
    peptide_table <- peptide_table[order(peptide_table$Log2_fc_delta, decreasing = FALSE),]
  }
  
  peptide_table <- within(peptide_table, e_coli_cumsum <- cumsum(Organism == "Escherichia coli (strain K12)"))
  peptide_table <- within(peptide_table, human_cumsum <- cumsum(Organism == "Homo sapiens"))
  peptide_table$FDP <- peptide_table$human_cumsum / (peptide_table$human_cumsum + peptide_table$e_coli_cumsum)
  
  pruned_peptide_table <- peptide_table[peptide_table$FDP <= fdp, ]
  
  fc_delta =  pruned_peptide_table$Log2_fc_delta %>% tail(n=1)
  e_coli_count =  pruned_peptide_table$e_coli_cumsum %>% tail(n=1)
  if(length(e_coli_count) == 0) e_coli_count = 0
  
  relevant_cols <- c("Organism", "Log2_fc", "Log2_fc_delta")
  
  return(list("fc_delta" = fc_delta, "true_positives" = e_coli_count,
              "peptide_table" = peptide_table))
}

all_conditions_iq <- c("A_", "B_", "C_",
                       "D_", "E_")
spike_in_levels <- c(1, 1.5, 2, 2.5, 3)
k <- 1
y <- 4
fc <- spike_in_levels[y] / spike_in_levels[k]


# With PIP
ions <- read.csv(iq_pip, sep = '\t')
ions[grepl("ECOLI", ions$Protein) & !grepl("HUMAN", ions$Mapped.Proteins), "Organism"] <- "Escherichia coli (strain K12)"
ions[grepl("HUMAN", ions$Protein) & !grepl("ECOLI", ions$Mapped.Proteins), "Organism"] <- "Homo sapiens"
ions <- subset(ions, Organism %in% 
                 c("Escherichia coli (strain K12)", 
                   "Homo sapiens"))

pip_results <- count_peps_at_fdp_iq(ions, conditions = c(all_conditions_iq[k], all_conditions_iq[y]), 
                                    fc = fc,
                                    fdp = 0.05, 
                                    order_by_fc = FALSE)

# w/o PIP
ions <- read.csv(iq_no_pip, sep = '\t')
ions[grepl("ECOLI", ions$Protein) & !grepl("HUMAN", ions$Mapped.Proteins), "Organism"] <- "Escherichia coli (strain K12)"
ions[grepl("HUMAN", ions$Protein) & !grepl("ECOLI", ions$Mapped.Proteins), "Organism"] <- "Homo sapiens"
ions <- subset(ions, Organism %in% 
                 c("Escherichia coli (strain K12)", 
                   "Homo sapiens"))
nopip_results <- count_peps_at_fdp_iq(ions, conditions = c(all_conditions_iq[k], all_conditions_iq[y]), 
                                      fc = fc,
                                      fdp = 0.05, 
                                      order_by_fc = FALSE)

# combine results
pip_results$peptide_table$OrgSetting <- paste0(pip_results$peptide_table$Organism, "_w_PIP")
nopip_results$peptide_table$OrgSetting <- paste0(nopip_results$peptide_table$Organism, "_no_PIP")

merged_tables <- rbind(pip_results$peptide_table, nopip_results$peptide_table)

# plot combined results
ymax <- 2000
text_height <- 70
bar_height <- 55
pep_text_height <- 185

ax_text_size <- 6
ax_title_size <- 8

y_max_zoom <- 225

full_plot <- ggplot(data = merged_tables, aes(fill = OrgSetting, x = Log2_fc)) +
  geom_histogram(alpha = 0.44, position = "identity", binwidth = 0.02, show.legend = F) +
  coord_cartesian(xlim = c(-2, 4), ylim = c(-2, ymax)) +
  scale_fill_manual(
    name = "Species",
    values = c("Homo sapiens_w_PIP" = "black", "Escherichia coli (strain K12)_w_PIP" = "#009970",
               "Homo sapiens_no_PIP" = "black", "Escherichia coli (strain K12)_no_PIP" = "#009970"),
    labels = c("H. sapiens", "E.coli", "H. sapiens", "E.coli")) +
  xlab("Log2 Fold Change") +
  ylab("Number of Peptides") +
  theme_bw() +
  theme(
        legend.position = "none",
        axis.text = element_text(size=ax_text_size),
        axis.title=element_text(size=ax_title_size))

ggsave(filename = paste0(r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\SubPlots_Revision\)", 
                         "IonQuant_Peptides_Full",
                         ".png"),
       plot = full_plot,
       units = "px",
       width = 1200,
       height = 1000)


bar_height <- 72
ymax <- 200
zoom_plot <- ggplot(data = merged_tables, aes(fill = OrgSetting, x = Log2_fc)) +
  geom_histogram(alpha = 0.44, position = "identity", binwidth = 0.02, show.legend = F) +
  coord_cartesian(xlim = c(-2, 4), ylim = c(-2, y_max_zoom)) +
  scale_fill_manual(
    name = "Species",
    values = c("Homo sapiens_w_PIP" = "black", "Escherichia coli (strain K12)_w_PIP" = "#009970",
               "Homo sapiens_no_PIP" = "black", "Escherichia coli (strain K12)_no_PIP" = "#009970"),
    labels = c("H. sapiens", "E.coli", "H. sapiens", "E.coli")) +
  xlab("Log2 Fold Change") +
  ylab("Number of Peptides") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=ax_text_size),
    axis.title=element_text(size=ax_title_size)) +
  geom_segment(y=0, yend=bar_height, x = log2(fc), xend = log2(fc), linetype="dashed", linewidth = 0.8) +
  geom_segment(y=0, yend=bar_height, x = log2(3.4), xend = log2(3.4), linetype="dashed", linewidth = 0.8, color = "red") 

ggsave(filename = paste0(r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\SubPlots_Revision\)", 
                         "IonQuant_Peptides_Zoom",
                         ".png"),
       plot = zoom_plot,
       units = "px",
       width = 1200,
       height = 1000)

# Specific Comparison FlashLFQ v1.0----

get_peptide_table <- function(path, condition_ids)
{
  peptides <- read.table(path, header = TRUE, sep = '\t')
  peptides <- peptides[!grepl("DECOY", peptides$Protein.Groups),]
  peptides <- subset(peptides, Organism %in% 
                       c("Escherichia coli (strain K12)", 
                         "Homo sapiens"))
  
  label_columns <- colnames(peptides)[c(1,3,4,5)]
  groupa <- get_col_names(peptides, condition_ids[1])
  groupb <-  get_col_names(peptides, condition_ids[2])
  groupa <- groupa[!sapply(groupa, function(x) grepl("Detection", x))]
  groupb <- groupb[!sapply(groupb, function(x) grepl("Detection", x))]
  
  peptides <- peptides[c(label_columns, groupa, groupb)]
  
  peptides[peptides==0] <- NA
  print(c(groupa, groupb))
  #peptides[c(groupa, groupb)] <- log2(peptides[c(groupa, groupb)])
  
  return(peptides)
}

count_peps_at_fdp <- function(peptide_table, conditions = c("A_", "E_"), fc = 3, fdp = 0.05, order_by_fc = F, sanity_check = F)
{
  cols <- colnames(peptide_table)[grepl(paste(conditions, collapse = "|"), colnames(peptide_table))]
  intensity_cols = cols[grepl("Intensity",cols)]
  i_cols_a <- intensity_cols[grepl(conditions[1], intensity_cols)]
  i_cols_e <- intensity_cols[grepl(conditions[2], intensity_cols)]
  
  # Mix them up such that the expected number of discoveries is zero
  if(sanity_check)
  {
    temp <- i_cols_a[c(3,4)]
    i_cols_a[c(3,4)] <- i_cols_e[c(3,4)]
    i_cols_e[c(3,4)] <- temp
  }
  
  peptide_table["A_median"] <- apply(peptide_table[i_cols_a], MARGIN = 1, function(x) mean(log2(x), na.rm=TRUE))
  peptide_table["E_median"] <- apply(peptide_table[i_cols_e], MARGIN = 1, function(x) mean(log2(x), na.rm=TRUE))
  peptide_table["Log2_fc"] <- peptide_table["E_median"] - peptide_table["A_median"]
  peptide_table <- peptide_table[!is.na(peptide_table$Log2_fc),]
  peptide_table["Log2_fc_delta"] <- abs(peptide_table$Log2_fc - log2(fc))
  
  if(order_by_fc)
  {
    peptide_table <- peptide_table[order(peptide_table$Log2_fc, decreasing = TRUE),]
  } else
  {
    peptide_table <- peptide_table[order(peptide_table$Log2_fc_delta, decreasing = FALSE),]
  }
  
  
  peptide_table <- within(peptide_table, e_coli_cumsum <- cumsum(Organism == "Escherichia coli (strain K12)"))
  peptide_table <- within(peptide_table, human_cumsum <- cumsum(Organism == "Homo sapiens"))
  peptide_table$FDP <- peptide_table$human_cumsum / (peptide_table$human_cumsum + peptide_table$e_coli_cumsum)
  
  pruned_peptide_table <- peptide_table[peptide_table$FDP <= fdp, ]
  fc_delta =  pruned_peptide_table$Log2_fc_delta %>% tail(n=1)
  e_coli_count =  pruned_peptide_table$e_coli_cumsum %>% tail(n=1)
  
  relevant_cols <- c("Organism", "Log2_fc", "Log2_fc_delta")
  
  return(list("fc_delta" = fc_delta, "true_positives" = e_coli_count,
              "peptide_table" = peptide_table))
}

get_fc_hist_specific <- function(results, comp, fc = 2.5, e_coli_count = 0, fc_delta = 0, title = "")
{
  if(length(e_coli_count) == 0)
  {
    e_coli_count <- 0
  }
  
  plot <- ggplot(data = results, aes(fill = Organism, x = Log2_fc)) +
    geom_histogram(alpha = 0.44, position = "identity", binwidth = 0.02, show.legend = F) +
    coord_cartesian(xlim = c(0, 2.5), ylim = c(-2, ymax)) +
    scale_fill_manual(
      name = "Species",
      values = c("Homo sapiens" = "black", "Escherichia coli (strain K12)" = "#009970"),
      labels = c("H. sapiens", "E.coli")) +
    xlab("Log2 Fold Change") +
    ylab("Number of Peptides") +
    theme_bw() +
    annotate("text", x = log2(fc), y = text_height, 
             label = paste0("E. coli peptides \nat 5% FDP: ", e_coli_count)) +
    annotate("text", x = -0.9, y = pep_text_height, 
             label = paste0("Total peptides \nquantified: ", nrow(results))) +
    geom_segment(y=0, yend=bar_height, x = log2(fc), xend = log2(fc), linetype="dashed", linewidth = 0.5)
  
  
  if(length(fc_delta) > 0 && fc_delta > 0)
  {
    plot <- plot + 
      geom_segment(y=0, yend=bar_height, x = log2(fc) - fc_delta, xend = log2(fc) - fc_delta, linetype="dashed", linewidth = 0.25) +
      geom_segment(y=0, yend=bar_height, x = log2(fc) + fc_delta, xend = log2(fc) + fc_delta, linetype="dashed", linewidth = 0.25)
  }
  
  title_size <- 16
  
  plot <- plot + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, size = title_size))
  
  return(plot)
}
all_conditions <- c("_A_", "_B_", "_C_",
                    "_D_", "_E_")
spike_in_levels <- c(1, 1.5, 2, 2.5, 3)
k <- 1
y <- 4
fc <- 2.5 / 1.0

old_flash_w_pep_path <- r"(D:\PXD003881_IonStar_SpikeIn\FlashLFQ_Current_Normed_NewPepQ\QuantifiedPeptides.tsv)"
ymax <- 200
text_height <- 185
bar_height <- 170

conditions <- c(all_conditions[k], all_conditions[y])
pip_table <- get_peptide_table(old_flash_w_pep_path, conditions)
fv1_w_pip <- count_peps_at_fdp(pip_table, conditions, fc, order_by_fc = FALSE, sanity_check = F)
fv1_w_pip_plot <- get_fc_hist_specific(fv1_w_pip$peptide_table, comparison, fc,
                                       e_coli_count = fv1_w_pip$true_positives, fc_delta = fv1_w_pip$fc_delta,
                                       title = "FlashLFQ v1.0 with PIP - 2.5x vs 1x")


old_flash_wo_pep_path <- r"(D:\PXD003881_IonStar_SpikeIn\FlashLFQ_Current_Normed_NewPepQ_NoMBR\QuantifiedPeptides.tsv)"

pip_table <- get_peptide_table(old_flash_wo_pep_path, conditions)
fv1_wo_pip <- count_peps_at_fdp(pip_table, conditions, fc, order_by_fc = FALSE, sanity_check = F)
fv1_wo_pip_plot <- get_fc_hist_specific(fv1_wo_pip$peptide_table, comparison, fc,
                                        e_coli_count = fv1_wo_pip$true_positives, fc_delta = fv1_wo_pip$fc_delta,
                                        title = "FlashLFQ v1.0 w/o PIP - 2.5x vs 1x")

# combine results
fv1_w_pip$peptide_table$OrgSetting <- paste0(fv1_w_pip$peptide_table$Organism, "_w_PIP")
fv1_wo_pip$peptide_table$OrgSetting <- paste0(fv1_wo_pip$peptide_table$Organism, "_no_PIP")

merged_tables_fv1 <- rbind(fv1_w_pip$peptide_table, fv1_wo_pip$peptide_table)



# plot combined results
ymax <- 2000
text_height <- 70
bar_height <- 55
pep_text_height <- 185

ax_text_size <- 6
ax_title_size <- 8

plot <- ggplot(data = merged_tables_fv1, aes(fill = OrgSetting, x = Log2_fc)) +
  geom_histogram(alpha = 0.44, position = "identity", binwidth = 0.02, show.legend = F) +
  coord_cartesian(xlim = c(-2, 4), ylim = c(-2, ymax)) +
  scale_fill_manual(
    name = "Species",
    values = c("Homo sapiens_w_PIP" = "black", "Escherichia coli (strain K12)_w_PIP" = "#009970",
               "Homo sapiens_no_PIP" = "black", "Escherichia coli (strain K12)_no_PIP" = "#009970"),
    labels = c("H. sapiens", "E.coli", "H. sapiens", "E.coli")) +
  xlab("Log2 Fold Change") +
  ylab("Number of Peptides") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=ax_text_size),
    axis.title=element_text(size=ax_title_size))

ggsave(filename = paste0(r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\SubPlots_Revision\)", 
                         "Flash_v1_Peptides_Full",
                         ".png"),
       plot = plot,
       units = "px",
       width = 1200,
       height = 1000)
plot


zoom_plot <- ggplot(data = merged_tables_fv1, aes(fill = OrgSetting, x = Log2_fc)) +
  geom_histogram(alpha = 0.44, position = "identity", binwidth = 0.02, show.legend = F) +
  coord_cartesian(xlim = c(-2, 4), ylim = c(-2, y_max_zoom)) +
  scale_fill_manual(
    name = "Species",
    values = c("Homo sapiens_w_PIP" = "black", "Escherichia coli (strain K12)_w_PIP" = "#009970",
               "Homo sapiens_no_PIP" = "black", "Escherichia coli (strain K12)_no_PIP" = "#009970"),
    labels = c("H. sapiens", "E.coli", "H. sapiens", "E.coli")) +
  xlab("Log2 Fold Change") +
  ylab("Number of Peptides") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=ax_text_size),
    axis.title=element_text(size=ax_title_size)) + 
  geom_segment(y=0, yend=172, x = log2(fc), xend = log2(fc), linetype="dashed", linewidth = 0.8)

ggsave(filename = paste0(r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\SubPlots_Revision\)", 
                         "Flash_v1_Peptides_zoom",
                         ".png"),
       plot = zoom_plot,
       units = "px",
       width = 1200,
       height = 1000)

# Specific Comparison FlashLFQ v2.0----

get_fc_hist_specific <- function(results, comp, fc = 2.5, e_coli_count = 0, fc_delta = 0, title = "")
{
  if(length(e_coli_count) == 0)
  {
    e_coli_count <- 0
  }
  
  plot <- ggplot(data = results, aes(fill = Organism, x = Log2_fc)) +
    geom_histogram(alpha = 0.44, position = "identity", binwidth = 0.02, show.legend = T) +
    coord_cartesian(xlim = c(0, 2.5), ylim = c(-2, ymax)) +
    scale_fill_manual(
      name = "Species",
      values = c("Homo sapiens" = "black", "Escherichia coli (strain K12)" = "#009970"),
      labels = c("E. coli", "H. sapiens")) +
    xlab("Log2 Fold Change") +
    ylab("Number of Peptides") +
    theme_bw()  +
    theme(panel.grid.major.y = element_line(color = "darkgray"),
          legend.position = "none",
          axis.text = element_text(size=12),
          axis.title=element_text(size=12),
          plot.title = element_text(size=12),
          strip.text.x = element_text(size = 12)) +
    annotate("text", x = 0.4, y = 165, 
             label = paste0("E. coli peptides \nat 5% FDP: ", e_coli_count), size = 4.5, hjust = 0) +
    annotate("text", x = 0.4, y = 140, 
             label = paste0("FC window width: ", format(fc_delta*2, digits=3)), size = 4.5, hjust = 0) +
    geom_segment(y=0, yend=bar_height, x = log2(fc), xend = log2(fc), linetype="dashed", linewidth = 0.5)
  
  
  if(length(fc_delta) > 0 && fc_delta > 0)
  {
    plot <- plot + 
      geom_segment(y=0, yend=bar_height, x = log2(fc) - fc_delta, xend = log2(fc) - fc_delta, linetype="dashed", linewidth = 0.25) +
      geom_segment(y=0, yend=bar_height, x = log2(fc) + fc_delta, xend = log2(fc) + fc_delta, linetype="dashed", linewidth = 0.25)
  }
  
  title_size <- 16
  
  
  plot <- plot + ggtitle(title) + theme(
    plot.title = element_text(hjust = 0.5, size = title_size),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 14))
  
  
  return(plot)
}

all_conditions <- c("_A_", "_B_", "_C_",
                    "_D_", "_E_")
spike_in_levels <- c(1, 1.5, 2, 2.5, 3)
k <- 1
y <- 4
fc <- 2.5 / 1.0

new_flash_w_pep_path <- r"(D:\PXD003881_IonStar_SpikeIn\FlashLFQ_7772_Normed_Donor1_Pip5\QuantifiedPeptides.tsv)"
ymax <- 200
text_height <- 185
bar_height <- 170


conditions <- c(all_conditions[k], all_conditions[y])
pip_table <- get_peptide_table(new_flash_w_pep_path, conditions)
fv2_w_pip <- count_peps_at_fdp(pip_table, conditions, fc, order_by_fc = FALSE, sanity_check = F)


# combine results
fv2_w_pip$peptide_table$OrgSetting <- paste0(fv2_w_pip$peptide_table$Organism, "_w_PIP")
fv1_wo_pip$peptide_table$OrgSetting <- paste0(fv1_wo_pip$peptide_table$Organism, "_no_PIP")

merged_tables_fv2 <- rbind(fv2_w_pip$peptide_table, fv1_wo_pip$peptide_table)

# plot combined results
ymax <- 2000
text_height <- 70
bar_height <- 55
pep_text_height <- 185

ax_text_size <- 6
ax_title_size <- 8

plot <- ggplot(data = merged_tables_fv2, aes(fill = OrgSetting, x = Log2_fc)) +
  geom_histogram(alpha = 0.44, position = "identity", binwidth = 0.02, show.legend = F) +
  coord_cartesian(xlim = c(-2, 4), ylim = c(-2, ymax)) +
  scale_fill_manual(
    name = "Species",
    values = c("Homo sapiens_w_PIP" = "black", "Escherichia coli (strain K12)_w_PIP" = "#009970",
               "Homo sapiens_no_PIP" = "black", "Escherichia coli (strain K12)_no_PIP" = "#009970"),
    labels = c("H. sapiens", "E.coli", "H. sapiens", "E.coli")) +
  xlab("Log2 Fold Change") +
  ylab("Number of Peptides") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=ax_text_size),
    axis.title=element_text(size=ax_title_size))

ggsave(filename = paste0(r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\SubPlots_Revision\)", 
                         "Flash_PipEcho_Peptides_Full",
                         ".png"),
       plot = plot,
       units = "px",
       width = 1200,
       height = 1000)
plot


zoom_plot <- ggplot(data = merged_tables_fv2, aes(fill = OrgSetting, x = Log2_fc)) +
  geom_histogram(alpha = 0.44, position = "identity", binwidth = 0.02, show.legend = F) +
  coord_cartesian(xlim = c(-2, 4), ylim = c(-2, y_max_zoom)) +
  scale_fill_manual(
    name = "Species",
    values = c("Homo sapiens_w_PIP" = "black", "Escherichia coli (strain K12)_w_PIP" = "#009970",
               "Homo sapiens_no_PIP" = "black", "Escherichia coli (strain K12)_no_PIP" = "#009970"),
    labels = c("H. sapiens", "E.coli", "H. sapiens", "E.coli")) +
  xlab("Log2 Fold Change") +
  ylab("Number of Peptides") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=ax_text_size),
    axis.title=element_text(size=ax_title_size)) + 
  geom_segment(y=0, yend=172, x = log2(fc), xend = log2(fc), linetype="dashed", linewidth = 0.8)

ggsave(filename = paste0(r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\SubPlots_Revision\)", 
                         "Flash_PipEcho_Peptides_zoom",
                         ".png"),
       plot = zoom_plot,
       units = "px",
       width = 1200,
       height = 1000)



#Specific Comparison MaxQuant ----
count_peps_at_fdp_mq <- function(peptide_table, conditions = c("A", "E"), fc = 3, fdp = 0.05, 
                                 protein_level = F)
{
  cols <- colnames(peptide_table)[grepl(paste(conditions, collapse = "|"), colnames(peptide_table))]
  if(protein_level)
  {
    intensity_cols = cols[grepl("LFQ.intensity",cols)]
  } else {
    intensity_cols = cols[grepl("Intensity",cols)]
  }
  
  i_cols_a <- intensity_cols[grepl(conditions[1], intensity_cols)]
  i_cols_e <- intensity_cols[grepl(conditions[2], intensity_cols)]
  
  peptide_table <- peptide_table[c("Organism", intensity_cols)]
  
  peptide_table[peptide_table == 0] <- NA

  peptide_table["A_median"] <- apply(peptide_table[i_cols_a], MARGIN = 1, function(x) mean(log2(x), na.rm=TRUE))
  peptide_table["E_median"] <- apply(peptide_table[i_cols_e], MARGIN = 1, function(x) mean(log2(x), na.rm=TRUE))
  peptide_table["Log2_fc"] <- peptide_table["E_median"] - peptide_table["A_median"]
  peptide_table <- peptide_table[!is.na(peptide_table$Log2_fc),]
  peptide_table["Log2_fc_delta"] <- abs(peptide_table$Log2_fc - log2(fc))

  peptide_table <- peptide_table[order(peptide_table$Log2_fc_delta, decreasing = FALSE),]
  
  peptide_table <- within(peptide_table, e_coli_cumsum <- cumsum(Organism == "Ecoli"))
  peptide_table <- within(peptide_table, human_cumsum <- cumsum(Organism == "Human"))
  
  peptide_table$FDP <- peptide_table$human_cumsum / (peptide_table$human_cumsum + peptide_table$e_coli_cumsum)
  
  last_rows <- which(peptide_table$FDP <= fdp, arr.ind = T)
  last_index = 1
  if(length(last_rows) > 0) 
  { 
    last_index = max(last_rows)
  }
  pruned_peptide_table <- peptide_table[1:last_index, ]
  
  fc_delta =  pruned_peptide_table$Log2_fc_delta %>% tail(n=1)
  e_coli_count =  pruned_peptide_table$e_coli_cumsum %>% tail(n=1)
  
  relevant_cols <- c("Organism", "Log2_fc", "Log2_fc_delta")
  
  return(list("fc_delta" = fc_delta, "true_positives" = e_coli_count,
              "peptide_table" = peptide_table[relevant_cols]))
}

mq_pep_paths = c(r"(D:\PXD003881_IonStar_SpikeIn\ReviewerResponse\MaxQuant_no_MBR_Ionstar\txt\peptides_Organism.txt)",
                 r"(D:\PXD003881_IonStar_SpikeIn\ReviewerResponse\MaxQuant_w_MBR_Ionstar\txt\peptides_Organism.txt)")
mq_conditions <- c("\\.A", "\\.B", "\\.C", "\\.D", "\\.E" )
spike_in_levels <- c(1, 1.5, 2, 2.5, 3)
k <- 1
y <- 4
fc <- 2.5 / 1.0
conditions <- c(mq_conditions[1], mq_conditions[4])

peptides <- read.table(mq_pep_paths[1], sep = '\t', header = T)
peptide_table<- subset(peptides, Organism %in% 
                         c("Ecoli", 
                           "Human"))
mq_no_mbr <- count_peps_at_fdp_mq(peptide_table, conditions, fc)

peptides <- read.table(mq_pep_paths[2], sep = '\t', header = T)
peptide_table<- subset(peptides, Organism %in% 
                         c("Ecoli", 
                           "Human"))
mq_w_mbr <- count_peps_at_fdp_mq(peptide_table, conditions, fc)



# combine results
mq_w_mbr$peptide_table$OrgSetting <- paste0(mq_w_mbr$peptide_table$Organism, "_w_PIP")
mq_no_mbr$peptide_table$OrgSetting <- paste0(mq_no_mbr$peptide_table$Organism, "_no_PIP")

merged_tables_mq <- rbind(mq_w_mbr$peptide_table, mq_no_mbr$peptide_table)

# plot combined results
ymax <- 2000
text_height <- 70
bar_height <- 55
pep_text_height <- 185

ax_text_size <- 6
ax_title_size <- 8

plot <- ggplot(data = merged_tables_mq, aes(fill = OrgSetting, x = Log2_fc)) +
  geom_histogram(alpha = 0.44, position = "identity", binwidth = 0.02, show.legend = F) +
  coord_cartesian(xlim = c(-2, 4), ylim = c(-2, ymax)) +
  scale_fill_manual(
    name = "Species",
    values = c("Human_w_PIP" = "black", "Ecoli_w_PIP" = "#009970",
               "Human_no_PIP" = "black", "Ecoli_no_PIP" = "#009970"),
    labels = c("H. sapiens", "E.coli", "H. sapiens", "E.coli")) +
  xlab("Log2 Fold Change") +
  ylab("Number of Peptides") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=ax_text_size),
    axis.title=element_text(size=ax_title_size))

ggsave(filename = paste0(r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\SubPlots_Revision\)", 
                         "MaxQuant_Peptides_Full",
                         ".png"),
       plot = plot,
       units = "px",
       width = 1200,
       height = 1000)
plot


zoom_plot <- ggplot(data = merged_tables_mq, aes(fill = OrgSetting, x = Log2_fc)) +
  geom_histogram(alpha = 0.44, position = "identity", binwidth = 0.02, show.legend = F) +
  coord_cartesian(xlim = c(-2, 4), ylim = c(-2, y_max_zoom)) +
  scale_fill_manual(
    name = "Species",
    values = c("Human_w_PIP" = "black", "Ecoli_w_PIP" = "#009970",
               "Human_no_PIP" = "black", "Ecoli_no_PIP" = "#009970"),
    labels = c("H. sapiens", "E.coli", "H. sapiens", "E.coli")) +
  xlab("Log2 Fold Change") +
  ylab("Number of Peptides") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=ax_text_size),
    axis.title=element_text(size=ax_title_size)) + 
  geom_segment(y=0, yend=95, x = log2(fc), xend = log2(fc), linetype="dashed", linewidth = 0.8)

ggsave(filename = paste0(r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\SubPlots_Revision\)", 
                         "MaxQuant_Peptides_zoom",
                         ".png"),
       plot = zoom_plot,
       units = "px",
       width = 1200,
       height = 1000)



# Protein Level Comparisons - IonQuant ----

# Paths
iq_pip <- "D:/PXD003881_IonStar_SpikeIn/IonQuant_Normed_Pip5/combined_protein.tsv"
iq_no_pip <- "D:/PXD003881_IonStar_SpikeIn/IonQuant_Normed_NoPip/combined_protein.tsv"


# Special functions
count_pros_at_fdp_iq <- function(peptide_table, conditions = c("A_", "E_"), fc = 3, fdp = 0.05)
{
  cols <- colnames(peptide_table)[grepl(paste(conditions, collapse = "|"), colnames(peptide_table))]
  intensity_cols = cols[grepl("Intensity",cols)]
  intensity_cols = intensity_cols[grepl("Max",intensity_cols)]
  i_cols_a <- intensity_cols[grepl(conditions[1], intensity_cols)]
  i_cols_e <- intensity_cols[grepl(conditions[2], intensity_cols)]
  peptide_table <- peptide_table[c("Organism", intensity_cols)]
  peptide_table[peptide_table == 0] <- NA
  
  peptide_table["A_median"] <- apply(peptide_table[i_cols_a], MARGIN = 1, function(x) mean(log2(x), na.rm=TRUE))
  peptide_table["E_median"] <- apply(peptide_table[i_cols_e], MARGIN = 1, function(x) mean(log2(x), na.rm=TRUE))
  peptide_table["Log2_fc"] <- peptide_table["E_median"] - peptide_table["A_median"]
  peptide_table <- peptide_table[!is.na(peptide_table$Log2_fc),]
  peptide_table["Log2_fc_delta"] <- abs(peptide_table$Log2_fc - log2(fc))
  
  peptide_table <- peptide_table[order(peptide_table$Log2_fc_delta, decreasing = FALSE),]
  
  peptide_table <- within(peptide_table, e_coli_cumsum <- cumsum(Organism == "Escherichia coli (strain K12)"))
  peptide_table <- within(peptide_table, human_cumsum <- cumsum(Organism == "Homo sapiens"))
  peptide_table$FDP <- peptide_table$human_cumsum / (peptide_table$human_cumsum + peptide_table$e_coli_cumsum)
  
  last_rows <- which(peptide_table$FDP <= fdp, arr.ind = T)
  last_index = 1
  if(length(last_rows) > 0) 
  { 
    last_index = max(last_rows)
  }
  pruned_peptide_table <- peptide_table[1:last_index, ]
  
  fc_delta =  pruned_peptide_table$Log2_fc_delta %>% tail(n=1)
  e_coli_count =  pruned_peptide_table$e_coli_cumsum %>% tail(n=1)
  if(length(e_coli_count) == 0) {
    e_coli_count = 0
  }
  
  relevant_cols <- c("Organism", "Log2_fc", "Log2_fc_delta", "FDP")
  
  return(list("fc_delta" = fc_delta, "true_positives" = e_coli_count,
              "peptide_table" = peptide_table[relevant_cols]))
}


all_conditions_iq <- c("A_", "B_", "C_",
                       "D_", "E_")
spike_in_levels <- c(1, 1.5, 2, 2.5, 3)
k <- 1
y <- 4
fc <- spike_in_levels[y] / spike_in_levels[k]


# With PIP
ions <- read.csv(iq_pip, sep = '\t')
ions[grepl("ECOLI", ions$Protein) & !grepl("HUMAN", ions$Mapped.Proteins), "Organism"] <- "Escherichia coli (strain K12)"
ions[grepl("HUMAN", ions$Protein) & !grepl("ECOLI", ions$Mapped.Proteins), "Organism"] <- "Homo sapiens"
ions <- subset(ions, Organism %in% 
                 c("Escherichia coli (strain K12)", 
                   "Homo sapiens"))

pip_results <- count_pros_at_fdp_iq(ions, conditions = c(all_conditions_iq[k], all_conditions_iq[y]), 
                                    fc = fc,
                                    fdp = 0.05)

# w/o PIP
ions <- read.csv(iq_no_pip, sep = '\t')
ions[grepl("ECOLI", ions$Protein) & !grepl("HUMAN", ions$Mapped.Proteins), "Organism"] <- "Escherichia coli (strain K12)"
ions[grepl("HUMAN", ions$Protein) & !grepl("ECOLI", ions$Mapped.Proteins), "Organism"] <- "Homo sapiens"
ions <- subset(ions, Organism %in% 
                 c("Escherichia coli (strain K12)", 
                   "Homo sapiens"))
nopip_results <- count_pros_at_fdp_iq(ions, conditions = c(all_conditions_iq[k], all_conditions_iq[y]), 
                                      fc = fc,
                                      fdp = 0.05)

# combine results
pip_results$peptide_table$OrgSetting <- paste0(pip_results$peptide_table$Organism, "_w_PIP")
nopip_results$peptide_table$OrgSetting <- paste0(nopip_results$peptide_table$Organism, "_no_PIP")

merged_tables <- rbind(pip_results$peptide_table, nopip_results$peptide_table)

# plot combined results
ymax <- 350
text_height <- 70
bar_height <- 55
pep_text_height <- 185

ax_text_size <- 6
ax_title_size <- 8

full_plot <- ggplot(data = merged_tables, aes(fill = OrgSetting, x = Log2_fc)) +
  geom_histogram(alpha = 0.44, position = "identity", binwidth = 0.02, show.legend = F) +
  coord_cartesian(xlim = c(-2, 4), ylim = c(-2, ymax)) +
  scale_fill_manual(
    name = "Species",
    values = c("Homo sapiens_w_PIP" = "black", "Escherichia coli (strain K12)_w_PIP" = "#009970",
               "Homo sapiens_no_PIP" = "black", "Escherichia coli (strain K12)_no_PIP" = "#009970"),
    labels = c("H. sapiens", "E.coli", "H. sapiens", "E.coli")) +
  xlab("Log2 Fold Change") +
  ylab("Number of Proteins") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=ax_text_size),
    axis.title=element_text(size=ax_title_size))
full_plot

ggsave(filename = paste0(r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\SubPlots_Revision\)", 
                         "IonQuant_Proteins_Full",
                         ".png"),
       plot = full_plot,
       units = "px",
       width = 1200,
       height = 1000)


bar_height <- 20
ymax <-75
zoom_plot <- ggplot(data = merged_tables, aes(fill = OrgSetting, x = Log2_fc)) +
  geom_histogram(alpha = 0.44, position = "identity", binwidth = 0.02, show.legend = F) +
  coord_cartesian(xlim = c(-1, 3), ylim = c(-2, ymax)) +
  scale_fill_manual(
    name = "Species",
    values = c("Homo sapiens_w_PIP" = "black", "Escherichia coli (strain K12)_w_PIP" = "#009970",
               "Homo sapiens_no_PIP" = "black", "Escherichia coli (strain K12)_no_PIP" = "#009970"),
    labels = c("H. sapiens", "E.coli", "H. sapiens", "E.coli")) +
  xlab("Log2 Fold Change") +
  ylab("Number of Proteins") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=ax_text_size),
    axis.title=element_text(size=ax_title_size)) +
  geom_segment(y=0, yend=bar_height, x = log2(fc), xend = log2(fc), linetype="dashed", linewidth = 0.8) +
  geom_segment(y=0, yend=bar_height, x = log2(3.5), xend = log2(3.5), linetype="dashed", linewidth = 0.8, color = "red") 

ggsave(filename = paste0(r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\SubPlots_Revision\)", 
                         "IonQuant_Proteins_Zoom",
                         ".png"),
       plot = zoom_plot,
       units = "px",
       width = 1200,
       height = 1000)

# Specific Comparison FlashLFQ v1.0----

get_peptide_table <- function(path, condition_ids)
{
  peptides <- read.table(path, header = TRUE, sep = '\t')
  peptides <- peptides[!grepl("DECOY", peptides$Protein.Groups),]
  peptides <- subset(peptides, Organism %in% 
                       c("Escherichia coli (strain K12)", 
                         "Homo sapiens"))
  
  label_columns <- colnames(peptides)[c(1,3,4,5)]
  groupa <- get_col_names(peptides, condition_ids[1])
  groupb <-  get_col_names(peptides, condition_ids[2])
  groupa <- groupa[!sapply(groupa, function(x) grepl("Detection", x))]
  groupb <- groupb[!sapply(groupb, function(x) grepl("Detection", x))]
  
  peptides <- peptides[c(label_columns, groupa, groupb)]
  
  peptides[peptides==0] <- NA
  print(c(groupa, groupb))
  #peptides[c(groupa, groupb)] <- log2(peptides[c(groupa, groupb)])
  
  return(peptides)
}

count_peps_at_fdp <- function(peptide_table, conditions = c("A_", "E_"), fc = 3, fdp = 0.05, order_by_fc = F, sanity_check = F)
{
  cols <- colnames(peptide_table)[grepl(paste(conditions, collapse = "|"), colnames(peptide_table))]
  intensity_cols = cols[grepl("Intensity",cols)]
  i_cols_a <- intensity_cols[grepl(conditions[1], intensity_cols)]
  i_cols_e <- intensity_cols[grepl(conditions[2], intensity_cols)]
  
  # Mix them up such that the expected number of discoveries is zero
  if(sanity_check)
  {
    temp <- i_cols_a[c(3,4)]
    i_cols_a[c(3,4)] <- i_cols_e[c(3,4)]
    i_cols_e[c(3,4)] <- temp
  }
  
  peptide_table["A_median"] <- apply(peptide_table[i_cols_a], MARGIN = 1, function(x) mean(log2(x), na.rm=TRUE))
  peptide_table["E_median"] <- apply(peptide_table[i_cols_e], MARGIN = 1, function(x) mean(log2(x), na.rm=TRUE))
  peptide_table["Log2_fc"] <- peptide_table["E_median"] - peptide_table["A_median"]
  peptide_table <- peptide_table[!is.na(peptide_table$Log2_fc),]
  peptide_table["Log2_fc_delta"] <- abs(peptide_table$Log2_fc - log2(fc))
  
  if(order_by_fc)
  {
    peptide_table <- peptide_table[order(peptide_table$Log2_fc, decreasing = TRUE),]
  } else
  {
    peptide_table <- peptide_table[order(peptide_table$Log2_fc_delta, decreasing = FALSE),]
  }
  
  
  peptide_table <- within(peptide_table, e_coli_cumsum <- cumsum(Organism == "Escherichia coli (strain K12)"))
  peptide_table <- within(peptide_table, human_cumsum <- cumsum(Organism == "Homo sapiens"))
  peptide_table$FDP <- peptide_table$human_cumsum / (peptide_table$human_cumsum + peptide_table$e_coli_cumsum)
  
  pruned_peptide_table <- peptide_table[peptide_table$FDP <= fdp, ]
  fc_delta =  pruned_peptide_table$Log2_fc_delta %>% tail(n=1)
  e_coli_count =  pruned_peptide_table$e_coli_cumsum %>% tail(n=1)
  
  relevant_cols <- c("Organism", "Log2_fc", "Log2_fc_delta")
  
  return(list("fc_delta" = fc_delta, "true_positives" = e_coli_count,
              "peptide_table" = peptide_table))
}

all_conditions <- c("_A_", "_B_", "_C_",
                    "_D_", "_E_")
spike_in_levels <- c(1, 1.5, 2, 2.5, 3)
k <- 1
y <- 4
fc <- 2.5 / 1.0

old_flash_w_pep_path <- r"(D:\PXD003881_IonStar_SpikeIn\FlashLFQ_Current_Normed_NewPepQ\QuantifiedProteins.tsv)"
ymax <- 200
text_height <- 185
bar_height <- 170

conditions <- c(all_conditions[k], all_conditions[y])
pip_table <- get_peptide_table(old_flash_w_pep_path, conditions)
fv1_w_pip <- count_peps_at_fdp(pip_table, conditions, fc, order_by_fc = FALSE, sanity_check = F)


old_flash_wo_pep_path <- r"(D:\PXD003881_IonStar_SpikeIn\FlashLFQ_Current_Normed_NewPepQ_NoMBR\QuantifiedProteins.tsv)"

pip_table <- get_peptide_table(old_flash_wo_pep_path, conditions)
fv1_wo_pip <- count_peps_at_fdp(pip_table, conditions, fc, order_by_fc = FALSE, sanity_check = F)


# combine results
fv1_w_pip$peptide_table$OrgSetting <- paste0(fv1_w_pip$peptide_table$Organism, "_w_PIP")
fv1_wo_pip$peptide_table$OrgSetting <- paste0(fv1_wo_pip$peptide_table$Organism, "_no_PIP")

merged_tables_fv1 <- rbind(fv1_w_pip$peptide_table, fv1_wo_pip$peptide_table)



# plot combined results
ymax <- 450
text_height <- 70
bar_height <- 55
pep_text_height <- 185

ax_text_size <- 6
ax_title_size <- 8

plot <- ggplot(data = merged_tables_fv1, aes(fill = OrgSetting, x = Log2_fc)) +
  geom_histogram(alpha = 0.44, position = "identity", binwidth = 0.02, show.legend = F) +
  coord_cartesian(xlim = c(-2, 4), ylim = c(-2, ymax)) +
  scale_fill_manual(
    name = "Species",
    values = c("Homo sapiens_w_PIP" = "black", "Escherichia coli (strain K12)_w_PIP" = "#009970",
               "Homo sapiens_no_PIP" = "black", "Escherichia coli (strain K12)_no_PIP" = "#009970"),
    labels = c("H. sapiens", "E.coli", "H. sapiens", "E.coli")) +
  xlab("Log2 Fold Change") +
  ylab("Number of Proteins") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=ax_text_size),
    axis.title=element_text(size=ax_title_size))

ggsave(filename = paste0(r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\SubPlots_Revision\)", 
                         "Flash_v1_Proteins_Full",
                         ".png"),
       plot = plot,
       units = "px",
       width = 1200,
       height = 1000)
plot


zoom_plot <- ggplot(data = merged_tables_fv1, aes(fill = OrgSetting, x = Log2_fc)) +
  geom_histogram(alpha = 0.44, position = "identity", binwidth = 0.02, show.legend = F) +
  coord_cartesian(xlim = c(-1, 3), ylim = c(-2, 75)) +
  scale_fill_manual(
    name = "Species",
    values = c("Homo sapiens_w_PIP" = "black", "Escherichia coli (strain K12)_w_PIP" = "#009970",
               "Homo sapiens_no_PIP" = "black", "Escherichia coli (strain K12)_no_PIP" = "#009970"),
    labels = c("H. sapiens", "E.coli", "H. sapiens", "E.coli")) +
  xlab("Log2 Fold Change") +
  ylab("Number of Proteins") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=ax_text_size),
    axis.title=element_text(size=ax_title_size)) + 
  geom_segment(y=0, yend=48, x = log2(fc), xend = log2(fc), linetype="dashed", linewidth = 0.8)

ggsave(filename = paste0(r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\SubPlots_Revision\)", 
                         "Flash_v1_Proteins_zoom",
                         ".png"),
       plot = zoom_plot,
       units = "px",
       width = 1200,
       height = 1000)

# Specific Comparison FlashLFQ v2.0----

all_conditions <- c("_A_", "_B_", "_C_",
                    "_D_", "_E_")
spike_in_levels <- c(1, 1.5, 2, 2.5, 3)
k <- 1
y <- 4
fc <- 2.5 / 1.0

new_flash_w_pep_path <- r"(D:\PXD003881_IonStar_SpikeIn\FlashLFQ_7772_Normed_Donor1_Pip5\QuantifiedProteins.tsv)"
ymax <- 200
text_height <- 185
bar_height <- 170

conditions <- c(all_conditions[k], all_conditions[y])
pip_table <- get_peptide_table(new_flash_w_pep_path, conditions)
fv2_w_pip <- count_peps_at_fdp(pip_table, conditions, fc, order_by_fc = FALSE, sanity_check = F)


# combine results
fv2_w_pip$peptide_table$OrgSetting <- paste0(fv2_w_pip$peptide_table$Organism, "_w_PIP")
fv1_wo_pip$peptide_table$OrgSetting <- paste0(fv1_wo_pip$peptide_table$Organism, "_no_PIP")

merged_tables_fv2 <- rbind(fv2_w_pip$peptide_table, fv1_wo_pip$peptide_table)

# plot combined results
ymax <- 450
text_height <- 70
bar_height <- 55
pep_text_height <- 185

ax_text_size <- 6
ax_title_size <- 8

plot <- ggplot(data = merged_tables_fv2, aes(fill = OrgSetting, x = Log2_fc)) +
  geom_histogram(alpha = 0.44, position = "identity", binwidth = 0.02, show.legend = F) +
  coord_cartesian(xlim = c(-2, 4), ylim = c(-2, ymax)) +
  scale_fill_manual(
    name = "Species",
    values = c("Homo sapiens_w_PIP" = "black", "Escherichia coli (strain K12)_w_PIP" = "#009970",
               "Homo sapiens_no_PIP" = "black", "Escherichia coli (strain K12)_no_PIP" = "#009970"),
    labels = c("H. sapiens", "E.coli", "H. sapiens", "E.coli")) +
  xlab("Log2 Fold Change") +
  ylab("Number of Proteins") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=ax_text_size),
    axis.title=element_text(size=ax_title_size))

ggsave(filename = paste0(r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\SubPlots_Revision\)", 
                         "Flash_PipEcho_Proteins_Full",
                         ".png"),
       plot = plot,
       units = "px",
       width = 1200,
       height = 1000)
plot


zoom_plot <- ggplot(data = merged_tables_fv2, aes(fill = OrgSetting, x = Log2_fc)) +
  geom_histogram(alpha = 0.44, position = "identity", binwidth = 0.02, show.legend = F) +
  coord_cartesian(xlim = c(-1, 3), ylim = c(-2, 75)) +
  scale_fill_manual(
    name = "Species",
    values = c("Homo sapiens_w_PIP" = "black", "Escherichia coli (strain K12)_w_PIP" = "#009970",
               "Homo sapiens_no_PIP" = "black", "Escherichia coli (strain K12)_no_PIP" = "#009970"),
    labels = c("H. sapiens", "E.coli", "H. sapiens", "E.coli")) +
  xlab("Log2 Fold Change") +
  ylab("Number of Proteins") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=ax_text_size),
    axis.title=element_text(size=ax_title_size)) + 
  geom_segment(y=0, yend=48, x = log2(fc), xend = log2(fc), linetype="dashed", linewidth = 0.8)

ggsave(filename = paste0(r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\SubPlots_Revision\)", 
                         "Flash_PipEcho_Proteins_zoom",
                         ".png"),
       plot = zoom_plot,
       units = "px",
       width = 1200,
       height = 1000)



#Specific Comparison MaxQuant ----
count_peps_at_fdp_mq <- function(peptide_table, conditions = c("A", "E"), fc = 3, fdp = 0.05, 
                                 protein_level = F)
{
  cols <- colnames(peptide_table)[grepl(paste(conditions, collapse = "|"), colnames(peptide_table))]
  if(protein_level)
  {
    intensity_cols = cols[grepl("LFQ.intensity",cols)]
  } else {
    intensity_cols = cols[grepl("Intensity",cols)]
  }
  
  i_cols_a <- intensity_cols[grepl(conditions[1], intensity_cols)]
  i_cols_e <- intensity_cols[grepl(conditions[2], intensity_cols)]
  
  peptide_table <- peptide_table[c("Organism", intensity_cols)]
  
  peptide_table[peptide_table == 0] <- NA
  
  peptide_table["A_median"] <- apply(peptide_table[i_cols_a], MARGIN = 1, function(x) mean(log2(x), na.rm=TRUE))
  peptide_table["E_median"] <- apply(peptide_table[i_cols_e], MARGIN = 1, function(x) mean(log2(x), na.rm=TRUE))
  peptide_table["Log2_fc"] <- peptide_table["E_median"] - peptide_table["A_median"]
  peptide_table <- peptide_table[!is.na(peptide_table$Log2_fc),]
  peptide_table["Log2_fc_delta"] <- abs(peptide_table$Log2_fc - log2(fc))
  
  peptide_table <- peptide_table[order(peptide_table$Log2_fc_delta, decreasing = FALSE),]
  
  peptide_table <- within(peptide_table, e_coli_cumsum <- cumsum(Organism == "Ecoli"))
  peptide_table <- within(peptide_table, human_cumsum <- cumsum(Organism == "Human"))
  
  peptide_table$FDP <- peptide_table$human_cumsum / (peptide_table$human_cumsum + peptide_table$e_coli_cumsum)
  
  last_rows <- which(peptide_table$FDP <= fdp, arr.ind = T)
  last_index = 1
  if(length(last_rows) > 0) 
  { 
    last_index = max(last_rows)
  }
  pruned_peptide_table <- peptide_table[1:last_index, ]
  
  fc_delta =  pruned_peptide_table$Log2_fc_delta %>% tail(n=1)
  e_coli_count =  pruned_peptide_table$e_coli_cumsum %>% tail(n=1)
  
  relevant_cols <- c("Organism", "Log2_fc", "Log2_fc_delta")
  
  return(list("fc_delta" = fc_delta, "true_positives" = e_coli_count,
              "peptide_table" = peptide_table[relevant_cols]))
}

mq_pep_paths = c(r"(D:\PXD003881_IonStar_SpikeIn\ReviewerResponse\MaxQuant_no_MBR_Ionstar\txt\proteinGroups_Organism.txt)",
                 r"(D:\PXD003881_IonStar_SpikeIn\ReviewerResponse\MaxQuant_w_MBR_Ionstar\txt\proteinGroups_Organism.txt)")

mq_conditions <- c("\\.A", "\\.B", "\\.C", "\\.D", "\\.E" )
spike_in_levels <- c(1, 1.5, 2, 2.5, 3)
k <- 1
y <- 4
fc <- 2.5 / 1.0
conditions <- c(mq_conditions[1], mq_conditions[4])

peptides <- read.table(mq_pep_paths[1], sep = '\t', header = T)
peptide_table<- subset(peptides, Organism %in% 
                         c("Ecoli", 
                           "Human"))
mq_no_mbr <- count_peps_at_fdp_mq(peptide_table, conditions, fc, protein_level = T)

peptides <- read.table(mq_pep_paths[2], sep = '\t', header = T)
peptide_table<- subset(peptides, Organism %in% 
                         c("Ecoli", 
                           "Human"))
mq_w_mbr <- count_peps_at_fdp_mq(peptide_table, conditions, fc, protein_level = T)



# combine results
mq_w_mbr$peptide_table$OrgSetting <- paste0(mq_w_mbr$peptide_table$Organism, "_w_PIP")
mq_no_mbr$peptide_table$OrgSetting <- paste0(mq_no_mbr$peptide_table$Organism, "_no_PIP")

merged_tables_mq <- rbind(mq_w_mbr$peptide_table, mq_no_mbr$peptide_table)

# plot combined results
ymax <- 450
text_height <- 70
bar_height <- 55
pep_text_height <- 185

ax_text_size <- 6
ax_title_size <- 8

plot <- ggplot(data = merged_tables_mq, aes(fill = OrgSetting, x = Log2_fc)) +
  geom_histogram(alpha = 0.44, position = "identity", binwidth = 0.02, show.legend = F) +
  coord_cartesian(xlim = c(-2, 4), ylim = c(-2, ymax)) +
  scale_fill_manual(
    name = "Species",
    values = c("Human_w_PIP" = "black", "Ecoli_w_PIP" = "#009970",
               "Human_no_PIP" = "black", "Ecoli_no_PIP" = "#009970"),
    labels = c("H. sapiens", "E.coli", "H. sapiens", "E.coli")) +
  xlab("Log2 Fold Change") +
  ylab("Number of Proteins") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=ax_text_size),
    axis.title=element_text(size=ax_title_size))

ggsave(filename = paste0(r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\SubPlots_Revision\)", 
                         "MaxQuant_Proteins_Full",
                         ".png"),
       plot = plot,
       units = "px",
       width = 1200,
       height = 1000)
plot


zoom_plot <- ggplot(data = merged_tables_mq, aes(fill = OrgSetting, x = Log2_fc)) +
  geom_histogram(alpha = 0.44, position = "identity", binwidth = 0.02, show.legend = F) +
  coord_cartesian(xlim = c(-1, 3), ylim = c(-2, 75)) +
  scale_fill_manual(
    name = "Species",
    values = c("Human_w_PIP" = "black", "Ecoli_w_PIP" = "#009970",
               "Human_no_PIP" = "black", "Ecoli_no_PIP" = "#009970"),
    labels = c("H. sapiens", "E.coli", "H. sapiens", "E.coli")) +
  xlab("Log2 Fold Change") +
  ylab("Number of Proteins") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=ax_text_size),
    axis.title=element_text(size=ax_title_size)) + 
  geom_segment(y=0, yend=36, x = log2(fc), xend = log2(fc), linetype="dashed", linewidth = 0.8)

ggsave(filename = paste0(r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\SubPlots_Revision\)", 
                         "MaxQuant_Proteins_zoom",
                         ".png"),
       plot = zoom_plot,
       units = "px",
       width = 1200,
       height = 1000)



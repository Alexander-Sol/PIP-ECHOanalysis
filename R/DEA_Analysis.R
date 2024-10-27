library(dplyr)
library(ggplot2)
library(limma)
library(tibble)
library(matrixStats)
library(scales)
library(patchwork)
library(cowplot)
library(ggpubr)

# Define functions for FlashLFQ output ----

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
  peptides[c(groupa, groupb)] <- log2(peptides[c(groupa, groupb)])
  
  return(peptides)
}

get_col_names <- function(table, condition_id)
{
  return(colnames(table)[grep(condition_id, colnames(table))])
}

get_full_peptide_table <- function(path)
{
  peptides <- read.table(path, header = TRUE, sep = '\t')
  peptides <- peptides[!grepl("DECOY", peptides$Protein.Groups),]
  peptides <- subset(peptides, Organism %in% 
                       c("Escherichia coli (strain K12)", 
                         "Homo sapiens"))
  
  peptides[peptides==0] <- NA
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

get_fc_hist <- function(results, comp)
{
  e_coli_count <- 0
  if(length(results$true_positives) > 0)
  {
    e_coli_count <- results$true_positives
  }
  
  plot <- ggplot(data = results$peptide_table, aes(fill = Organism, x = Log2_fc)) +
    geom_histogram(alpha = 0.2, position = "identity", binwidth = 0.02, show.legend = F) +
    coord_cartesian(xlim = c(-0.5, 1.75), ylim = c(-10, 350)) +
    scale_fill_manual(
      name = "Species",
      values = c("Homo sapiens" = "black", "Escherichia coli (strain K12)" = "red"),
      labels = c("H. sapiens", "E.coli")) +
    xlab("Log2 Fold Change") +
    ylab("Number of Peptides") +
    theme_bw() +
    annotate("text", x = 0.8, y = 280, 
             label = paste0("E. coli peptides \nat 5% FDP: ", e_coli_count))
  
  title_size <- 16
  if (grepl("_A_", comp))
  {
    if(grepl("_B_", comp)) 
      plot <- plot + ggtitle("1.5x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    if(grepl("_C_", comp)) 
      plot <- plot + ggtitle("2x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    if(grepl("_D_", comp)) 
      plot <- plot + ggtitle("2.5x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    if(grepl("_E_", comp)) 
      plot <- plot + ggtitle("3x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
  }
  
  return(plot)
}

get_fc_hist_pro <- function(results, comp)
{
  e_coli_count <- 0
  if(length(results$true_positives) > 0)
  {
    e_coli_count <- results$true_positives
  }
  
  plot <- ggplot(data = results$peptide_table, aes(fill = Organism, x = Log2_fc)) +
    geom_histogram(alpha = 0.2, position = "identity", binwidth = 0.02, show.legend = F) +
    coord_cartesian(xlim = c(-0.5, 1.75), ylim = c(-5, 100)) +
    scale_fill_manual(
      name = "Species",
      values = c("Homo sapiens" = "black", "Escherichia coli (strain K12)" = "red"),
      labels = c("H. sapiens", "E.coli")) +
    xlab("Log2 Fold Change") +
    ylab("Number of Peptides") +
    theme_bw() +
    annotate("text", x = 0.8, y = 80, 
             label = paste0("E. coli proteins \nat 5% FDP: ", e_coli_count))
  
  title_size <- 16
  if (grepl("_A_", comp))
  {
    if(grepl("_B_", comp)) 
      plot <- plot + ggtitle("1.5x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    if(grepl("_C_", comp)) 
      plot <- plot + ggtitle("2x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    if(grepl("_D_", comp)) 
      plot <- plot + ggtitle("2.5x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    if(grepl("_E_", comp)) 
      plot <- plot + ggtitle("3x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
  }
  
  return(plot)
}

multi_comp_analysis <- function(path, protein_level = F, order_by_fc = F, sanity_check = F)
{
  pep_table <- get_full_peptide_table(path)
  
  all_conditions <- c("_A_", "_B_", "_C_",
                      "_D_", "_E_")
  spike_in_levels <- c(1, 1.5, 2, 2.5, 3)
  deltas <- c()
  tps <- c()
  comparison <- c()
  log_fcs <- c()
  combined_pep_table <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(combined_pep_table) <- c("Organism", "Log2_fc", "Log2_fc_delta", "Condition")
  hist_plots <- vector("list", 10)
  iter = 1
  
  for(k in 1:4)
  {
    for(y in (k+1):5)
    {
      
      fc <- (spike_in_levels[y]/ spike_in_levels[k])
      results <- count_peps_at_fdp(pep_table, conditions = c(all_conditions[k], all_conditions[y]), 
                                   fc = fc,
                                   fdp = 0.05,
                                   order_by_fc = order_by_fc,
                                   sanity_check = sanity_check)
      comp = paste0(all_conditions[k], "vs.", all_conditions[y])
      if(length(results$fc_delta) == 0)
      {
        deltas <- c(deltas, NA)
        tps <- c(tps, NA)
      }
      else
      {
        deltas <- c(deltas, results$fc_delta)
        tps <- c(tps, results$true_positives)
      }
      
      comparison <- c(comparison, comp)
      log_fcs <- c(log_fcs, fc)
      
      #Now we're going to merge all the results in the peptide table
      results$peptide_table$Condition <- comp
      combined_pep_table <- rbind(combined_pep_table, results$peptide_table)
      
      if(protein_level)
      {
        hist_plots[[iter]] <- get_fc_hist_pro(results, comp)
      } else
      {
        hist_plots[[iter]] <- get_fc_hist(results, comp)
      }
      
      iter <- iter + 1
    }
  }
  
  
  tps
  comparison
  
  total_peps <- ggplot() + annotate("text", x = 0, y = 0, label = 
                                      paste0("E. coli peptides at 5% FDP\nacross all comparisons: ",
                                             sum(tps, na.rm=T))) + theme_void()
  
  leg_plot <- ggplot(data = results$peptide_table, aes(fill = Organism, x = Log2_fc)) +
    geom_histogram(alpha = 0.2, position = "identity", binwidth = 0.02, show.legend = T) +
    scale_fill_manual(
      name = "Species",
      values = c("Homo sapiens" = "black", "Escherichia coli (strain K12)" = "red"),
      labels = c( "E.coli", "H. sapiens")) 
  
  legend <- get_legend(leg_plot)
  as_ggplot(legend)
  
  r_title_size <- 5
  r1 <- ggplot() + annotate("text", x = 0, y = 0, angle = -90, label = "1x", size = r_title_size) + theme_void()
  r2 <- ggplot() + annotate("text", x = 0, y = 0, angle = -90, label = "1.5x", size = r_title_size) + theme_void()
  r3 <- ggplot() + annotate("text", x = 0, y = 0, angle = -90, label = "2x", size = r_title_size) + theme_void()
  r4 <- ggplot() + annotate("text", x = 0, y = 0, angle = -90, label = "2.5x", size = r_title_size) + theme_void()
  
  gc()
  
  multi_plot <- plot_grid( 
    hist_plots[[1]], hist_plots[[2]], hist_plots[[3]],hist_plots[[4]], r1,
    NULL, hist_plots[[5]], hist_plots[[6]], hist_plots[[7]], r2,
    as_ggplot(legend), NULL, hist_plots[[8]], hist_plots[[9]], r3,
    total_peps, NULL, NULL, hist_plots[[10]], r4,
    ncol = 5,
    rel_widths = c(4,4,4,4,1))
  
  return(list("plot" = multi_plot, "tp_count" = sum(tps, na.rm=T)))
}

# MaxQuant Functions ----
count_peps_at_fdp_mq <- function(peptide_table, conditions = c("A", "E"), fc = 3, fdp = 0.05, order_by_fc = F, sanity_check = F,
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

get_fc_hist_mq <- function(results, comp)
{
  e_coli_count <- 0
  if(length(results$true_positives) > 0)
  {
    e_coli_count <- results$true_positives
  }
  
  plot <- ggplot(data = results$peptide_table, aes(fill = Organism, x = Log2_fc)) +
    geom_histogram(alpha = 0.2, position = "identity", binwidth = 0.02, show.legend = F) +
    coord_cartesian(xlim = c(-0.5, 1.75), ylim = c(-10, 350)) +
    scale_fill_manual(
      name = "Species",
      values = c("Human" = "black", "Ecoli" = "red"),
      labels = c("H. sapiens", "E.coli")) +
    xlab("Log2 Fold Change") +
    ylab("Number of Peptides") +
    theme_bw() +
    annotate("text", x = 0.8, y = 280, 
             label = paste0("E. coli peptides \nat 5% FDP: ", e_coli_count))
  
  title_size <- 16
  if (grepl("A", comp))
  {
    if(grepl("B", comp)) 
      plot <- plot + ggtitle("1.5x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    if(grepl("C", comp)) 
      plot <- plot + ggtitle("2x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    if(grepl("D", comp)) 
      plot <- plot + ggtitle("2.5x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    if(grepl("E", comp)) 
      plot <- plot + ggtitle("3x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
  }
  
  return(plot)
}

get_fc_hist_pro_mq <- function(results, comp)
{
  e_coli_count <- 0
  if(length(results$true_positives) > 0)
  {
    e_coli_count <- results$true_positives
  }
  
  plot <- ggplot(data = results$peptide_table, aes(fill = Organism, x = Log2_fc)) +
    geom_histogram(alpha = 0.2, position = "identity", binwidth = 0.02, show.legend = F) +
    coord_cartesian(xlim = c(-0.5, 1.75), ylim = c(-5, 100)) +
    scale_fill_manual(
      name = "Species",
      values = c("Human" = "black", "Ecoli" = "red"),
      labels = c("H. sapiens", "E.coli")) +
    xlab("Log2 Fold Change") +
    ylab("Number of Peptides") +
    theme_bw() +
    annotate("text", x = 0.8, y = 80, 
             label = paste0("E. coli proteins \nat 5% FDP: ", e_coli_count))
  
  title_size <- 16
  if (grepl("Av", comp))
  {
    if(grepl("Bv", comp)) 
      plot <- plot + ggtitle("1.5x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    if(grepl("Cv", comp)) 
      plot <- plot + ggtitle("2x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    if(grepl("Dv", comp)) 
      plot <- plot + ggtitle("2.5x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    if(grepl("Ev", comp)) 
      plot <- plot + ggtitle("3x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
  }
  
  return(plot)
}

multi_comp_analysis_mq <- function(path, protein_level = F, order_by_fc = F, sanity_check = F)
{
  peptides <- read.table(path, sep = '\t', header = T)
  pep_table<- subset(peptides, Organism %in% 
                       c("Ecoli", 
                         "Human"))
  
  all_conditions <- c("A", "B", "C", "D", "E")
  spike_in_levels <- c(1, 1.5, 2, 2.5, 3)
  deltas <- c()
  tps <- c()
  comparison <- c()
  log_fcs <- c()
  combined_pep_table <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(combined_pep_table) <- c("Organism", "Log2_fc", "Log2_fc_delta", "Condition")
  hist_plots <- vector("list", 10)
  iter = 1
  
  for(k in 1:4)
  {
    for(y in (k+1):5)
    {
      
      fc <- (spike_in_levels[y]/ spike_in_levels[k])
      results <- count_peps_at_fdp_mq(pep_table, conditions = c(all_conditions[k], all_conditions[y]), 
                                      fc = fc,
                                      fdp = 0.05,
                                      order_by_fc = order_by_fc,
                                      sanity_check = sanity_check,
                                      protein_level = protein_level)
      comp = paste0(all_conditions[k], "vs.", all_conditions[y])
      if(length(results$fc_delta) == 0)
      {
        deltas <- c(deltas, NA)
        tps <- c(tps, NA)
      }
      else
      {
        deltas <- c(deltas, results$fc_delta)
        tps <- c(tps, results$true_positives)
      }
      
      comparison <- c(comparison, comp)
      log_fcs <- c(log_fcs, fc)
      
      #Now we're going to merge all the results in the peptide table
      results$peptide_table$Condition <- comp
      combined_pep_table <- rbind(combined_pep_table, results$peptide_table)
      
      if(protein_level)
      {
        hist_plots[[iter]] <- get_fc_hist_pro_mq(results, comp)
      } else
      {
        hist_plots[[iter]] <- get_fc_hist_mq(results, comp)
      }
      
      iter <- iter + 1
    }
  }
  
  
  tps
  comparison
  
  total_peps <- ggplot() + annotate("text", x = 0, y = 0, label = 
                                      paste0("E. coli peptides at 5% FDP\nacross all comparisons: ",
                                             sum(tps, na.rm=T))) + theme_void()
  
  leg_plot <- ggplot(data = results$peptide_table, aes(fill = Organism, x = Log2_fc)) +
    geom_histogram(alpha = 0.2, position = "identity", binwidth = 0.02, show.legend = T) +
    scale_fill_manual(
      name = "Species",
      values = c("Human" = "black", "Ecoli" = "red"),
      labels = c( "E.coli", "H. sapiens")) 
  
  legend <- get_legend(leg_plot)
  as_ggplot(legend)
  
  r_title_size <- 5
  r1 <- ggplot() + annotate("text", x = 0, y = 0, angle = -90, label = "1x", size = r_title_size) + theme_void()
  r2 <- ggplot() + annotate("text", x = 0, y = 0, angle = -90, label = "1.5x", size = r_title_size) + theme_void()
  r3 <- ggplot() + annotate("text", x = 0, y = 0, angle = -90, label = "2x", size = r_title_size) + theme_void()
  r4 <- ggplot() + annotate("text", x = 0, y = 0, angle = -90, label = "2.5x", size = r_title_size) + theme_void()
  
  gc()
  
  multi_plot <- plot_grid( 
    hist_plots[[1]], hist_plots[[2]], hist_plots[[3]],hist_plots[[4]], r1,
    NULL, hist_plots[[5]], hist_plots[[6]], hist_plots[[7]], r2,
    as_ggplot(legend), NULL, hist_plots[[8]], hist_plots[[9]], r3,
    total_peps, NULL, NULL, hist_plots[[10]], r4,
    ncol = 5,
    rel_widths = c(4,4,4,4,1))
  
  return(list("plot" = multi_plot, "tp_count" = sum(tps, na.rm=T)))
}

# IonQuant functions ----

count_peps_at_fdp_iq <- function(peptide_table, conditions = c("A_", "E_"), fc = 3, fdp = 0.05, order_by_fc = F)
{
  cols <- colnames(peptide_table)[grepl(paste(conditions, collapse = "|"), colnames(peptide_table))]
  cols <- colnames(peptide_table)[grepl(paste(conditions, collapse = "|"), colnames(peptide_table))]
  intensity_cols = cols[grepl("Intensity",cols)]
  i_cols_a <- intensity_cols[grepl(conditions[1], intensity_cols)]
  i_cols_e <- intensity_cols[grepl(conditions[2], intensity_cols)]
  
  peptide_table$quant_count <- apply(peptide_table, 1, function(x) length(which(x %in% c("MBR", "MS/MS"))))
  peptide_table <- peptide_table %>% 
    arrange(desc(quant_count)) %>%
    group_by(Modified.Sequence) %>%
    slice(1)
  
  
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

get_fc_hist_iq <- function(results, comp)
{
  e_coli_count <- 0
  if(length(results$true_positives) > 0)
  {
    e_coli_count <- results$true_positives
  }
  
  plot <- ggplot(data = results$peptide_table, aes(fill = Organism, x = Log2_fc)) +
    geom_histogram(alpha = 0.2, position = "identity", binwidth = 0.02, show.legend = F) +
    coord_cartesian(xlim = c(-1, 4), ylim = c(-10, 200)) +
    scale_fill_manual(
      name = "Species",
      values = c("Homo sapiens" = "black", "Escherichia coli (strain K12)" = "red"),
      labels = c("H. sapiens", "E.coli")) +
    xlab("Log2 Fold Change") +
    ylab("Number of Peptides") +
    theme_bw() +
    annotate("text", x = 2, y = 150, 
             label = paste0("E. coli peptides \nat 5% FDP: ", e_coli_count))
  
  title_size <- 16
  if (grepl("A_", comp))
  {
    if(grepl("B_", comp)) 
      plot <- plot + ggtitle("1.5x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    if(grepl("C_", comp)) 
      plot <- plot + ggtitle("2x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    if(grepl("D_", comp)) 
      plot <- plot + ggtitle("2.5x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    if(grepl("E_", comp)) 
      plot <- plot + ggtitle("3x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
  }
  
  return(plot)
}

multi_comp_iq_pep <- function(ion_path, order_by_fc = F)
{
  ions <- read.csv(ion_path, sep = '\t')
  ions[grepl("ECOLI", ions$Protein) & !grepl("HUMAN", ions$Mapped.Proteins), "Organism"] <- "Escherichia coli (strain K12)"
  ions[grepl("HUMAN", ions$Protein) & !grepl("ECOLI", ions$Mapped.Proteins), "Organism"] <- "Homo sapiens"
  ions <- subset(ions, Organism %in% 
                   c("Escherichia coli (strain K12)", 
                     "Homo sapiens"))
  
  all_conditions_iq <- c("A_", "B_", "C_",
                         "D_", "E_")
  spike_in_levels <- c(1, 1.5, 2, 2.5, 3)
  deltas <- c()
  tps <- c()
  comparison <- c()
  log_fcs <- c()
  combined_pep_table <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(combined_pep_table) <- c("Organism", "Log2_fc", "Log2_fc_delta", "Condition")
  hist_plots <- vector("list", 10)
  iter = 1
  
  for(k in 1:4)
  {
    for(y in (k+1):5)
    {
      
      fc <- (spike_in_levels[y]/ spike_in_levels[k])
      results <- count_peps_at_fdp_iq(ions, conditions = c(all_conditions_iq[k], all_conditions_iq[y]), 
                                      fc = fc,
                                      fdp = 0.05, 
                                      order_by_fc = order_by_fc)
      comp = paste0(all_conditions_iq[k], "vs.", all_conditions_iq[y])
      if(length(results$fc_delta) == 0)
      {
        deltas <- c(deltas, NA)
        tps <- c(tps, NA)
      }
      else
      {
        deltas <- c(deltas, results$fc_delta)
        tps <- c(tps, results$true_positives)
      }
      
      comparison <- c(comparison, comp)
      log_fcs <- c(log_fcs, fc)
      
      #Now we're going to merge all the results in the peptide table
      results$peptide_table$Condition <- comp
      combined_pep_table <- rbind(combined_pep_table, results$peptide_table)
      
      hist_plots[[iter]] <- get_fc_hist_iq(results, comp)
      iter <- iter + 1
    }
  }
  
  total_peps <- ggplot() + annotate("text", x = 0, y = 0, label = 
                                      paste0("E. coli peptides at 5% FDP\nacross all comparisons: ",
                                             sum(tps, na.rm=T))) + theme_void()
  
  leg_plot <- ggplot(data = results$peptide_table, aes(fill = Organism, x = Log2_fc)) +
    geom_histogram(alpha = 0.2, position = "identity", binwidth = 0.02, show.legend = T) +
    scale_fill_manual(
      name = "Species",
      values = c("Homo sapiens" = "black", "Escherichia coli (strain K12)" = "red"),
      labels = c( "E.coli", "H. sapiens")) 
  
  legend <- get_legend(leg_plot)
  as_ggplot(legend)
  
  r_title_size <- 5
  r1 <- ggplot() + annotate("text", x = 0, y = 0, angle = -90, label = "1x", size = r_title_size) + theme_void()
  r2 <- ggplot() + annotate("text", x = 0, y = 0, angle = -90, label = "1.5x", size = r_title_size) + theme_void()
  r3 <- ggplot() + annotate("text", x = 0, y = 0, angle = -90, label = "2x", size = r_title_size) + theme_void()
  r4 <- ggplot() + annotate("text", x = 0, y = 0, angle = -90, label = "2.5x", size = r_title_size) + theme_void()
  
  gc()
  
  multi_plot <- plot_grid( 
    hist_plots[[1]], hist_plots[[2]], hist_plots[[3]],hist_plots[[4]], r1,
    NULL, hist_plots[[5]], hist_plots[[6]], hist_plots[[7]], r2,
    as_ggplot(legend), NULL, hist_plots[[8]], hist_plots[[9]], r3,
    total_peps, NULL, NULL, hist_plots[[10]], r4,
    ncol = 5,
    rel_widths = c(4,4,4,4,1))
  
  return(list("plot" = multi_plot, "tp_count" = sum(tps, na.rm=T)))
}

# Counting PIP IDs ----
# FlashLFQ was ran multiple times, and the different output were written to the same
# output directory, specified here
base_path = "D:/..." 
donors = c("Donor02", "Donor05","Donor1")
pips = c("_Pip1", "_Pip2p5", "_Pip5", "_Pip100")

base_out_dir <- r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\Peptide_Level_10_25_24\)"
dir.create(base_out_dir, showWarnings = FALSE)

# PIP-ECHO + FlashLFQ - Count MBR Results by Analysis Settings ----

file_path = "/QuantifiedPeptides.tsv"

all_peptide_paths = c()
pip_settings = c()
mbr_counts = c()
msms_counts = c()
for(i in 1:3)
{
  for(j in 1:4)
  {
    peptide_path = paste0(base_path, donors[i], pips[j], file_path)
    pip_settings <- c(pip_settings, paste0(donors[i], pips[j]))
    pep_res <- get_full_peptide_table(peptide_path)
    mbr_counts <- c(mbr_counts, sum(pep_res == "MBR", na.rm = TRUE))
    msms_counts <- c(msms_counts, sum(pep_res == "MSMS", na.rm = TRUE))
  }
}

# IonQuant ----

file_path = "/combined_ion.tsv"
base_path = "D:/PXD003881_IonStar_SpikeIn/IonQuant_Normed"

mbr_counts_iq = c()
msms_counts_iq = c()
pip_settings_iq = c()

for(j in 1:4)
{
  ion_path = paste0(base_path, pips[j], file_path)
  pip_settings_iq <- c(pip_settings_iq, pips[j])
  ions <- read.csv(ion_path, sep = '\t')
  ions$MBR_Count <- apply(ions, 1, function(x) length(which(x=="MBR")))
  ions <- ions %>% 
    arrange(desc(MBR_Count)) %>%
    group_by(Modified.Sequence) %>%
    slice(1)
  mbr_counts_iq <- c(mbr_counts_iq, sum(ions == "MBR", na.rm = TRUE))
  
  msms_counts_iq <- c(msms_counts_iq, sum(ions == "MS/MS", na.rm = TRUE))
}

mbr_counts_iq

# FlashLFQ v1.0 ----
old_peptide_path <- r"(D:\PXD003881_IonStar_SpikeIn\FlashLFQ_Current_Normed_MM106_wMBR\QuantifiedPeptides.tsv)"
peptide_path = old_peptide_path
pep_res <- get_full_peptide_table(peptide_path)
mbr_counts_old <- sum(pep_res == "MBR", na.rm = TRUE)
msms_counts_old <- sum(pep_res == "MSMS", na.rm = TRUE)

# MaxQuant ----
mq_path <- r"(D:\PXD003881_IonStar_SpikeIn\txt_MQ_ionQuant_spikeIn\evidence.txt)"
mq_res <- read.csv(mq_path, sep = '\t')
mbr_counts_mq = sum(mq_res$Type %in% c("MULTI-MATCH", "MULTI-MATCH-MSMS"), na.rm = T)
msms_counts_mq = sum(mq_res$Type %in% c("MULTI-MSMS", "MULTI-SECPEP"), na.rm = T)


# Plotting PIP Events ----

mbr_combined <- c(mbr_counts_old, mbr_counts_mq, mbr_counts_iq, mbr_counts)
msms_combined <- c(msms_counts_old, msms_counts_mq, msms_counts_iq, msms_counts)

pip_list <- list(pip_condition = c("FlashLFQ\nv1.0", "MaxQuant", paste0("IonQuant", pip_settings_iq), pip_settings),
                 pip_counts = mbr_combined,
                 msms_count = msms_combined)
pip_table <- as.data.frame(pip_list)
pip_table$pip_condition <- factor(pip_table$pip_condition,
                                  levels = c("FlashLFQ\nv1.0", "MaxQuant", 
                                             paste0("IonQuant", pip_settings_iq), pip_settings))
pip_table$DonorFilter <- factor(c("FlashLFQ\nv1.0", "MaxQuant",  rep("IonQuant", 4),
                                  rep("Donor FDR Cutoff = 0.2%", 4),
                                  rep("Donor FDR Cutoff = 0.5%", 4),
                                  rep("Donor FDR Cutoff = 1%", 4)),
                                levels = c("FlashLFQ\nv1.0", "MaxQuant", "IonQuant",
                                           "Donor FDR Cutoff = 0.2%",
                                           "Donor FDR Cutoff = 0.5%",
                                           "Donor FDR Cutoff = 1%"))

pip_table$PipFilter <- factor(c(100, 100, 1, 2.5, 5, 100, 1, 2.5, 5, 100, 1, 2.5, 5, 100, 1, 2.5, 5, 100))

# Here, select the data that will be plotted
# Each PIP-ECHO output FDR corresponds to a different donor peptide FDR

mask <- (pip_table$PipFilter == 1 & pip_table$DonorFilter == "Donor FDR Cutoff = 0.2%") |
  (pip_table$PipFilter == 2.5 & pip_table$DonorFilter == "Donor FDR Cutoff = 0.5%") |
  (pip_table$PipFilter == 5 & pip_table$DonorFilter == "Donor FDR Cutoff = 1%") |
  (pip_table$PipFilter == 100 & pip_table$DonorFilter == "Donor FDR Cutoff = 1%") |
  pip_table$DonorFilter == "FlashLFQ\nv1.0" |
  pip_table$DonorFilter == "IonQuant" |
  pip_table$DonorFilter ==  "MaxQuant"

pip_table_subset <- pip_table[mask,]

# rename and re-order 
levels(pip_table_subset$DonorFilter)[levels(pip_table_subset$DonorFilter) %in% 
                                       c("Donor FDR Cutoff = 0.2%", 
                                         "Donor FDR Cutoff = 0.5%",
                                         "Donor FDR Cutoff = 1%")] <- "FlashLFQ +\nPIP-ECHO"

pip_table_subset$DonorFilter <- factor(pip_table_subset$DonorFilter,
                                       levels = c("FlashLFQ +\nPIP-ECHO",
                                                  "IonQuant",
                                                  "MaxQuant",
                                                  "FlashLFQ\nv1.0"))

pip_plot <- ggplot(data= pip_table_subset, aes(x = PipFilter, y = pip_counts)) +
  geom_bar(stat="identity", position=position_dodge(), fill =  alpha("#0072B2", .75)) +
  facet_grid(~DonorFilter, space = "free", scales = "free") +
  ylab("Number of PIP Events") +
  xlab("PIP FDR Cutoff (%)") +
  scale_y_continuous(labels = comma) +
  theme_bw()  +
  theme(axis.text = element_text(size=12),
        axis.title=element_text(size=14))
pip_plot

ggsave(filename = paste0(base_out_dir, "PIP_Plot_sub", ".png"),
       plot = pip_plot,
       units = "px",
       width = 2700,
       height = 2100)


# Peptide Level Analysis - FlashFLQ ----

base_path = "D:/PXD003881_IonStar_SpikeIn/FlashLFQ_7772_Normed_"
donors = c("Donor02", "Donor05","Donor1")
pips = c("_Pip1", "_Pip2p5", "_Pip5", "_Pip100")
file_path = "/QuantifiedPeptides.tsv"

tp_counts <- c()
pip_settings <- c()
for(i in 1:3)
{
  for(j in 1:4)
  {
    peptide_path = paste0(base_path, donors[i], pips[j], file_path)
    multi_res <- multi_comp_analysis(peptide_path, order_by_fc = F)
    ggsave(filename = paste0(base_out_dir, donors[i], pips[j], ".png"),
           plot = multi_res$plot,
           units = "px",
           width = 3161,
           height = 2343)
    
    pip_settings <- c(pip_settings, paste0(donors[i], pips[j]))
    tp_counts <- c(tp_counts, multi_res$tp_count)
  }
}

tp_counts

# Peptide level FlashLFQ v1.0 Analysis ----

flash_v1_base <- r"(D:\PXD003881_IonStar_SpikeIn\FlashLFQ_Current_Normed_NewPepQ)"
pips = c("_NoMBR", "")
file_path = "/QuantifiedPeptides.tsv"
pip_settings_v1 = c()
tp_counts_v1 = c()

for(i in 1:2)
{
  pip_settings_v1 <- c(pip_settings_v1, pips[i])
  old_flash_path = paste0(flash_v1_base, pips[i], file_path)
  multi_res_v1_pep <- multi_comp_analysis(old_flash_path, order_by_fc = F)
  ggsave(filename = paste0(base_out_dir, "FlashLFQ_v1", pips[i], ".png"),
         plot = multi_res_v1_pep$plot,
         units = "px",
         width = 3161,
         height = 2343)
  tp_counts_v1 <- c(tp_counts_v1, multi_res_v1_pep$tp_count)
}


# Peptide level IonQuant Analysis ----

base_path_iq = "D:/PXD003881_IonStar_SpikeIn/IonQuant_Normed"
pips = c("_NoPip", "_Pip1", "_Pip2p5", "_Pip5", "_Pip100")
file_path_iq = "/combined_ion.tsv"

tp_counts_iq <- c()
pip_settings_iq <- c()

for(j in 1:5)
{
  ion_path = paste0(base_path_iq, pips[j], file_path_iq)
  
  multi_res <- multi_comp_iq_pep(ion_path)
  ggsave(filename = paste0(base_out_dir, "IonQuant_", pips[j], ".png"),
         plot = multi_res$plot,
         units = "px",
         width = 3161,
         height = 2343)
  
  pip_settings_iq <- c(pip_settings_iq, pips[j])
  tp_counts_iq <- c(tp_counts_iq, multi_res$tp_count)
}

tp_counts_iq

# Peptide level MaxQuant Analysis ----

mq_pep_paths <- c(r"(D:\PXD003881_IonStar_SpikeIn\combined_no_MBR_IonStar\txt\peptides_with_organisms.txt)",
                  r"(D:\PXD003881_IonStar_SpikeIn\txt_MQ_ionQuant_spikeIn\peptides_with_organisms.txt)")

mq_settings <- c("No PIP", "100")
tp_counts_mq <- c()
for (i in 1:2)
{
  multi_res_mq_pep <- multi_comp_analysis_mq(mq_pep_paths[i])
  ggsave(filename = paste0(base_out_dir, mq_settings[i], ".png"),
         plot = multi_res_mq_pep$plot,
         units = "px",
         width = 3161,
         height = 2343)
  tp_counts_mq <- c(tp_counts_mq, multi_res_mq_pep$tp_count)
}

# Plotting Peptide Level summary statistics ----

tp_counts_combined <- c(tp_counts_v1, tp_counts_mq, tp_counts_iq, tp_counts)

count_list <- list(pip_condition = c(paste0("FlashLFQ v1", pip_settings_v1),
                                     paste0("MaxQuant", mq_settings),
                                     paste0("IonQuant", pip_settings_iq),
                                     pip_settings),
                   counts = tp_counts_combined)
count_table_pep <- as.data.frame(count_list)
count_table_pep$pip_condition <- factor(count_table_pep$pip_condition,
                                        levels = c(paste0("FlashLFQ v1", pip_settings_v1),
                                                   paste0("MaxQuant", mq_settings), 
                                                   paste0("IonQuant", pip_settings_iq),
                                                   pip_settings))
count_table_pep$DonorFilter <- factor(c(rep("FlashLFQ v1", 2), "MaxQuant", "MaxQuant",
                                        rep("IonQuant", 5),
                                        rep("Donor FDR Cutoff = 0.2%", 4),
                                        rep("Donor FDR Cutoff = 0.5%", 4),
                                        rep("Donor FDR Cutoff = 1%", 4)),
                                      levels = c("FlashLFQ v1", "MaxQuant", "IonQuant",
                                                 "Donor FDR Cutoff = 0.2%",
                                                 "Donor FDR Cutoff = 0.5%",
                                                 "Donor FDR Cutoff = 1%"))

count_table_pep$PipFilter <- factor(c("No PIP", 100,
                                      "No PIP", 100,
                                      "No PIP", 1, 2.5, 5, 100,
                                      1, 2.5, 5, 100,
                                      1, 2.5, 5, 100,
                                      1, 2.5, 5, 100),
                                    levels = c("No PIP",1, 2.5, 5, 100))


# Select the condition we want to include in the main figure
# Each PIP-ECHO output FDR corresponds to a different donor peptide FDR
mask <- (count_table_pep$PipFilter == 1 & count_table_pep$DonorFilter == "Donor FDR Cutoff = 0.2%") |
  (count_table_pep$PipFilter == 2.5 & count_table_pep$DonorFilter == "Donor FDR Cutoff = 0.5%") |
  (count_table_pep$PipFilter == 5 & count_table_pep$DonorFilter == "Donor FDR Cutoff = 1%") |
  (count_table_pep$PipFilter == 100 & count_table_pep$DonorFilter == "Donor FDR Cutoff = 1%") |
  count_table_pep$DonorFilter == "FlashLFQ v1" |
  count_table_pep$DonorFilter == "IonQuant" |
  count_table_pep$DonorFilter ==  "MaxQuant"

count_table_pep_subset <- count_table_pep[mask,]

# rename and re-order 
levels(count_table_pep_subset$DonorFilter)[levels(count_table_pep_subset$DonorFilter) == "FlashLFQ v1"] <- "FlashLFQ\nv1.0"
levels(count_table_pep_subset$DonorFilter)[levels(count_table_pep_subset$DonorFilter) %in% 
                                             c("Donor FDR Cutoff = 0.2%", 
                                               "Donor FDR Cutoff = 0.5%",
                                               "Donor FDR Cutoff = 1%")] <- "FlashLFQ +\nPIP-ECHO"

count_table_pep_subset$DonorFilter <- factor(count_table_pep_subset$DonorFilter,
                                             levels = c("FlashLFQ +\nPIP-ECHO",
                                                        "IonQuant",
                                                        "MaxQuant",
                                                        "FlashLFQ\nv1.0"))
flash_no_pip <- count_table_pep_subset[count_table_pep_subset$DonorFilter == "FlashLFQ\nv1.0" 
                                       & count_table_pep_subset$PipFilter == "No PIP",]
flash_no_pip$DonorFilter <- "FlashLFQ +\nPIP-ECHO"
count_table_pep_subset <- rbind(count_table_pep_subset, flash_no_pip)

dea_plot_pep_sub <- ggplot(data = count_table_pep_subset, aes(x = PipFilter, y = counts)) +
  geom_bar(stat="identity", position=position_dodge(), fill = alpha("darkred", .75)) +
  facet_grid(~DonorFilter, space = "free", scales = "free") +
  ylab("Differentially abundant E.coli peptides at 5% FDP") +
  xlab("PIP FDR Cutoff (%)") +
  ggtitle("Number of E. coli Peptides with Correct Differential Abundance") +
  scale_y_continuous(labels = comma) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title=element_text(size=14))
dea_plot_pep_sub


ggsave(filename = paste0(base_out_dir, "DEA_Plot_sub", ".png"),
       plot = dea_plot_pep_sub,
       units = "px",
       width = 2700,
       height = 2100)

# Protein Level Analysis ----

base_path = "D:/PXD003881_IonStar_SpikeIn/FlashLFQ_7772_Normed_"
donors = c("Donor02", "Donor05","Donor1")
pips = c("_Pip1", "_Pip2p5", "_Pip5", "_Pip100")
file_path = "/QuantifiedProteins.tsv"
base_out_dir <- r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\Protein_Level_10_25_24\)"
dir.create(base_out_dir, showWarnings = F)

tp_counts <- c()
pip_settings <- c()
for(i in 1:3)
{
  for(j in 1:4)
  {
    peptide_path = paste0(base_path, donors[i], pips[j], file_path)
    multi_res <- multi_comp_analysis(peptide_path, protein_level = T, order_by_fc = F)
    ggsave(filename = paste0(base_out_dir, donors[i], pips[j], ".png"),
           plot = multi_res$plot,
           units = "px",
           width = 3161,
           height = 2343)
    
    pip_settings <- c(pip_settings, paste0(donors[i], pips[j]))
    tp_counts <- c(tp_counts, multi_res$tp_count)
  }
}

tp_counts

# FlashLFQ v1.0 Protein level ----

flash_v1_base <- r"(D:\PXD003881_IonStar_SpikeIn\FlashLFQ_Current_Normed_NewPepQ)"
pips = c("_NoMBR", "")
file_path = "/QuantifiedProteins.tsv"
pip_settings_v1 = c()
tp_counts_v1 = c()

for(i in 1:2)
{
  pip_settings_v1 <- c(pip_settings_v1, pips[i])
  old_flash_path = paste0(flash_v1_base, pips[i], file_path)
  multi_res_v1_pep <- multi_comp_analysis(old_flash_path, protein_level = T)
  ggsave(filename = paste0(base_out_dir, "FlashLFQ_v1", pips[i], ".png"),
         plot = multi_res_v1_pep$plot,
         units = "px",
         width = 3161,
         height = 2343)
  tp_counts_v1 <- c(tp_counts_v1, multi_res_v1_pep$tp_count)
}
tp_counts_v1


# IonQuant Protein level analysis ----

count_pros_at_fdp_iq <- function(peptide_table, conditions = c("A_", "E_"), fc = 3, fdp = 0.05)
{
  cols <- colnames(peptide_table)[grepl(paste(conditions, collapse = "|"), colnames(peptide_table))]
  intensity_cols = cols[grepl("Intensity",cols)]
  intensity_cols = intensity_cols[grepl("Max",intensity_cols)]
  i_cols_a <- intensity_cols[grepl(conditions[1], intensity_cols)]
  i_cols_e <- intensity_cols[grepl(conditions[2], intensity_cols)]
  
  # peptide_table["A_median"] <- apply(peptide_table[i_cols_a], MARGIN = 1, function(x) log2(mean(x,na.rm=TRUE)))
  # peptide_table["E_median"] <- apply(peptide_table[i_cols_e], MARGIN = 1, function(x) log2(mean(x,na.rm=TRUE)))
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
  
  relevant_cols <- c("Organism", "Log2_fc", "Log2_fc_delta", "FDP")
  
  return(list("fc_delta" = fc_delta, "true_positives" = e_coli_count,
              "peptide_table" = peptide_table[relevant_cols]))
}

get_fc_hist_iq_pro <- function(results, comp)
{
  e_coli_count <- 0
  if(length(results$true_positives) > 0)
  {
    e_coli_count <- results$true_positives
  }
  
  plot <- ggplot(data = results$peptide_table, aes(fill = Organism, x = Log2_fc)) +
    geom_histogram(alpha = 0.2, position = "identity", binwidth = 0.02, show.legend = F) +
    coord_cartesian(xlim = c(-1, 6), ylim = c(-2, 50)) +
    scale_fill_manual(
      name = "Species",
      values = c("Homo sapiens" = "black", "Escherichia coli (strain K12)" = "red"),
      labels = c("H. sapiens", "E.coli")) +
    xlab("Log2 Fold Change") +
    ylab("Number of Proteins") +
    theme_bw() +
    annotate("text", x = 4, y = 25, 
             label = paste0("E. coli proteins \nat 5% FDP: ", e_coli_count))
  
  title_size <- 16
  if (grepl("A_", comp))
  {
    if(grepl("B_", comp)) 
      plot <- plot + ggtitle("1.5x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    if(grepl("C_", comp)) 
      plot <- plot + ggtitle("2x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    if(grepl("D_", comp)) 
      plot <- plot + ggtitle("2.5x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    if(grepl("E_", comp)) 
      plot <- plot + ggtitle("3x") + theme(plot.title = element_text(hjust = 0.5, size = title_size))
  }
  
  return(plot)
}

multi_comp_iq <- function(ion_path, protein_level = T)
{
  ions <- read.csv(ion_path, sep = '\t')
  ions[grepl("ECOLI", ions$Protein), "Organism"] <- "Escherichia coli (strain K12)"
  ions[grepl("HUMAN", ions$Protein), "Organism"] <- "Homo sapiens"
  ions <- subset(ions, Organism %in% 
                   c("Escherichia coli (strain K12)", 
                     "Homo sapiens"))
  
  all_conditions_iq <- c("A_", "B_", "C_",
                         "D_", "E_")
  spike_in_levels <- c(1, 1.5, 2, 2.5, 3)
  deltas <- c()
  tps <- c()
  comparison <- c()
  log_fcs <- c()
  combined_pep_table <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(combined_pep_table) <- c("Organism", "Log2_fc", "Log2_fc_delta", "Condition")
  hist_plots <- vector("list", 10)
  iter = 1
  
  for(k in 1:4)
  {
    for(y in (k+1):5)
    {
      
      fc <- (spike_in_levels[y]/ spike_in_levels[k])
      results <- count_pros_at_fdp_iq(ions, conditions = c(all_conditions_iq[k], all_conditions_iq[y]), 
                                      fc = fc,
                                      fdp = 0.05)
      comp = paste0(all_conditions_iq[k], "vs.", all_conditions_iq[y])
      if(length(results$fc_delta) == 0)
      {
        deltas <- c(deltas, NA)
        tps <- c(tps, NA)
      }
      else
      {
        deltas <- c(deltas, results$fc_delta)
        tps <- c(tps, results$true_positives)
      }
      
      comparison <- c(comparison, comp)
      log_fcs <- c(log_fcs, fc)
      
      #Now we're going to merge all the results in the peptide table
      results$peptide_table$Condition <- comp
      combined_pep_table <- rbind(combined_pep_table, results$peptide_table)
      
      hist_plots[[iter]] <- get_fc_hist_iq_pro(results, comp)
      iter <- iter + 1
    }
  }
  
  total_peps <- ggplot() + annotate("text", x = 0, y = 0, label = 
                                      paste0("E. coli proteins at 5% FDP\nacross all comparisons: ",
                                             sum(tps, na.rm=T))) + theme_void()
  
  leg_plot <- ggplot(data = results$peptide_table, aes(fill = Organism, x = Log2_fc)) +
    geom_histogram(alpha = 0.2, position = "identity", binwidth = 0.02, show.legend = T) +
    scale_fill_manual(
      name = "Species",
      values = c("Homo sapiens" = "black", "Escherichia coli (strain K12)" = "red"),
      labels = c( "E.coli", "H. sapiens")) 
  
  legend <- get_legend(leg_plot)
  as_ggplot(legend)
  
  r_title_size <- 5
  r1 <- ggplot() + annotate("text", x = 0, y = 0, angle = -90, label = "1x", size = r_title_size) + theme_void()
  r2 <- ggplot() + annotate("text", x = 0, y = 0, angle = -90, label = "1.5x", size = r_title_size) + theme_void()
  r3 <- ggplot() + annotate("text", x = 0, y = 0, angle = -90, label = "2x", size = r_title_size) + theme_void()
  r4 <- ggplot() + annotate("text", x = 0, y = 0, angle = -90, label = "2.5x", size = r_title_size) + theme_void()
  
  gc()
  
  multi_plot <- plot_grid( 
    hist_plots[[1]], hist_plots[[2]], hist_plots[[3]],hist_plots[[4]], r1,
    NULL, hist_plots[[5]], hist_plots[[6]], hist_plots[[7]], r2,
    as_ggplot(legend), NULL, hist_plots[[8]], hist_plots[[9]], r3,
    total_peps, NULL, NULL, hist_plots[[10]], r4,
    ncol = 5,
    rel_widths = c(4,4,4,4,1))
  
  return(list("plot" = multi_plot, "tp_count" = sum(tps, na.rm=T)))
}


base_path_iq = "D:/PXD003881_IonStar_SpikeIn/IonQuant_Normed"
pips = c("_NoPip", "_Pip1", "_Pip2p5", "_Pip5", "_Pip100")
file_path_iq = "/combined_protein.tsv"

tp_counts_iq <- c()
pip_settings_iq <- c()

for(j in 1:5)
{
  ion_path = paste0(base_path_iq, pips[j], file_path_iq)
  
  multi_res <- multi_comp_iq(ion_path, protein_level = T)
  ggsave(filename = paste0(base_out_dir, "IonQuant_", pips[j], ".png"),
         plot = multi_res$plot,
         units = "px",
         width = 3161,
         height = 2343)
  
  pip_settings_iq <- c(pip_settings_iq, pips[j])
  tp_counts_iq <- c(tp_counts_iq, multi_res$tp_count)
}

# MaxQuant Analysis ----

mq_pep_paths <- c(r"(D:\PXD003881_IonStar_SpikeIn\combined_no_MBR_IonStar\txt\proteinGroups_with_organisms.txt)",
                  r"(D:\PXD003881_IonStar_SpikeIn\txt_MQ_ionQuant_spikeIn\proteinGroups_with_organisms.txt)")

mq_settings <- c("No PIP", "100")
tp_counts_mq <- c()
for (i in 1:2)
{
  multi_res_mq_pep <- multi_comp_analysis_mq(mq_pep_paths[i], protein_level = T)
  ggsave(filename = paste0(base_out_dir, mq_settings[i], ".png"),
         plot = multi_res_mq_pep$plot,
         units = "px",
         width = 3161,
         height = 2343)
  tp_counts_mq <- c(tp_counts_mq, multi_res_mq_pep$tp_count)
}

# Protein summary statistics ----

tp_counts_combined <- c(tp_counts_v1, tp_counts_mq, tp_counts_iq, tp_counts)

count_list <- list(pip_condition = c(paste0("FlashLFQ v1", pip_settings_v1),
                                     paste0("MaxQuant", mq_settings), 
                                     paste0("IonQuant", pip_settings_iq),
                                     pip_settings),
                   counts = tp_counts_combined)
count_table <- as.data.frame(count_list)
count_table$pip_condition <-  factor(count_table$pip_condition,
                                     levels = c(paste0("FlashLFQ v1", pip_settings_v1),
                                                paste0("MaxQuant", mq_settings), 
                                                paste0("IonQuant", pip_settings_iq),
                                                pip_settings))
count_table$DonorFilter <- factor(c(rep("FlashLFQ v1", 2), "MaxQuant", "MaxQuant",
                                    rep("IonQuant", 5),
                                    rep("Donor FDR Cutoff = 0.2%", 4),
                                    rep("Donor FDR Cutoff = 0.5%", 4),
                                    rep("Donor FDR Cutoff = 1%", 4)),
                                  levels = c("FlashLFQ v1", "MaxQuant", "IonQuant",
                                             "Donor FDR Cutoff = 0.2%",
                                             "Donor FDR Cutoff = 0.5%",
                                             "Donor FDR Cutoff = 1%"))

count_table$PipFilter <- factor(c("No PIP", 100,
                                  "No PIP", 100,
                                  "No PIP", 1, 2.5, 5, 100,
                                  1, 2.5, 5, 100,
                                  1, 2.5, 5, 100,
                                  1, 2.5, 5, 100),
                                levels = c("No PIP",1, 2.5, 5, 100))


# Select the condition we want to include in the main figure
# Each PIP-ECHO output FDR corresponds to a different donor peptide FDR
mask <- (count_table$PipFilter == 1 & count_table$DonorFilter == "Donor FDR Cutoff = 0.2%") |
  (count_table$PipFilter == 2.5 & count_table$DonorFilter == "Donor FDR Cutoff = 0.5%") |
  (count_table$PipFilter == 5 & count_table$DonorFilter == "Donor FDR Cutoff = 1%") |
  (count_table$PipFilter == 100 & count_table$DonorFilter == "Donor FDR Cutoff = 1%") |
  count_table$DonorFilter == "FlashLFQ v1" |
  count_table$DonorFilter == "IonQuant" |
  count_table$DonorFilter ==  "MaxQuant"

count_table_subset <- count_table[mask,]

# rename and re-order 
levels(count_table_subset$DonorFilter)[levels(count_table_subset$DonorFilter) == "FlashLFQ v1"] <- "FlashLFQ\nv1.0"
levels(count_table_subset$DonorFilter)[levels(count_table_subset$DonorFilter) %in% 
                                         c("Donor FDR Cutoff = 0.2%", 
                                           "Donor FDR Cutoff = 0.5%",
                                           "Donor FDR Cutoff = 1%")] <- "FlashLFQ +\nPIP-ECHO"

count_table_subset$DonorFilter <- factor(count_table_subset$DonorFilter,
                                         levels = c("FlashLFQ +\nPIP-ECHO",
                                                    "IonQuant",
                                                    "MaxQuant",
                                                    "FlashLFQ\nv1.0"))
flash_no_pip <- count_table_subset[count_table_subset$DonorFilter == "FlashLFQ\nv1.0" 
                                   & count_table_subset$PipFilter == "No PIP",]
flash_no_pip$DonorFilter <- "FlashLFQ +\nPIP-ECHO"
count_table_subset <- rbind(count_table_subset, flash_no_pip)

dea_plot_pro_sub <- ggplot(data = count_table_subset, aes(x = PipFilter, y = counts)) +
  geom_bar(stat="identity", position=position_dodge(), fill = alpha("darkred", .75)) +
  facet_grid(~DonorFilter, space = "free", scales = "free") +
  ylab("Differentially abundant E.coli proteins at 5% FDP") +
  xlab("PIP FDR Cutoff (%)") +
  ggtitle("Number of E. coli Proteins with Correct Differential Abundance") +
  scale_y_continuous(labels = comma) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title=element_text(size=14))
dea_plot_pro_sub


ggsave(filename = paste0(base_out_dir, "DEA_Plot_pro_sub", ".png"),
       plot = dea_plot_pro_sub,
       units = "px",
       width = 2700,
       height = 2100)




count_table_sub <- count_table[which(!(count_table$DonorFilter %in% c("Donor FDR Cutoff = 0.2%", "Donor FDR Cutoff = 5%"))), ]
levels(count_table_sub$DonorFilter)[levels(count_table_sub$DonorFilter) == "FlashLFQ v1"] <- "FlashLFQ\nv1.0"
levels(count_table_sub$DonorFilter)[levels(count_table_sub$DonorFilter) == "Donor FDR Cutoff = 1%"] <- "FlashLFQ\nv2.0"

dea_plot_sub <- ggplot(data = count_table_sub, aes(x = PipFilter, y = counts)) +
  geom_bar(stat="identity", position=position_dodge(), fill = alpha("darkred", .75)) +
  facet_grid(~DonorFilter, space = "free", scales = "free") +
  ylab("Differentially abundant E.coli proteins at 5% FDP") +
  xlab("PIP FDR Cutoff (%)") +
  ggtitle("Number of E. coli Proteins with Correct Differential Abundance") +
  scale_y_continuous(labels = comma) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title=element_text(size=14))


ggsave(filename = paste0(base_out_dir, "DEA_Plot_sub", ".png"),
       plot = dea_plot_sub,
       units = "px",
       width = 2700,
       height = 2100)

# All Plotting ---- 

base_out_dir <- r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\SubPlots_10_25_24\)"
dir.create(base_out_dir, showWarnings = FALSE)

levels(count_table_subset$DonorFilter) <- c("FlashLFQ +\nPIP-ECHO","IonQuant","Max-\nQuant","FlashLFQ\nv1.0")

levels(count_table_pep_subset$DonorFilter) <- c("FlashLFQ +\nPIP-ECHO","IonQuant","Max-\nQuant","FlashLFQ\nv1.0")

ax_text_size = 15
ax_title_size = 20
plot_title_size = 20
strip_text_size = 13

dea_plot_pro_sub <- ggplot(data = count_table_subset, aes(x = PipFilter, y = counts)) +
  geom_bar(stat="identity", position=position_dodge(), fill = alpha("#009E73", .75)) +
  facet_grid(~DonorFilter, space = "free", scales = "free") +
  ylab("E.coli proteins at 5% FDP") +
  xlab("PIP FDR Cutoff (%)") +
  ggtitle("Number of Proteins with Correct Differential Abundance") +
  scale_y_continuous(labels = comma) +
  theme_bw() +
  theme(panel.grid.major.y = element_line(color = "darkgray"),
        legend.position = "none",
        axis.text = element_text(size=ax_text_size),
        axis.title=element_text(size=ax_title_size),
        plot.title = element_text(size=plot_title_size),
        strip.text.x = element_text(size = strip_text_size))


ggsave(filename = paste0(base_out_dir, "DEA_Plot_pro_sub", ".png"),
       plot = dea_plot_pro_sub,
       units = "px",
       width = 2700,
       height = 2100)


dea_plot_pep_sub <- ggplot(data = count_table_pep_subset, aes(x = PipFilter, y = counts)) +
  geom_bar(stat="identity", position=position_dodge(), fill = alpha("#009E73", .75)) +
  facet_grid(~DonorFilter, space = "free", scales = "free") +
  ylab("E.coli peptides at 5% FDP") +
  xlab("PIP FDR Cutoff (%)") +
  ggtitle("Number of Peptides with Correct Differential Abundance") +
  scale_y_continuous(labels = comma) +
  theme_bw() +
  theme(panel.grid.major.y = element_line(color = "darkgray"),
        legend.position = "none",
        axis.text = element_text(size=ax_text_size),
        axis.title=element_text(size=ax_title_size),
        plot.title = element_text(size=plot_title_size),
        strip.text.x = element_text(size = strip_text_size))


ggsave(filename = paste0(base_out_dir, "DEA_Plot_Pep_sub", ".png"),
       plot = dea_plot_pep_sub,
       units = "px",
       width = 2700,
       height = 2100)

ax_text_size = 15
ax_title_size = 20
plot_title_size = 24
strip_text_size = 13

levels(pip_table_subset$DonorFilter) <- c("FlashLFQ +\nPIP-ECHO","IonQuant","Max-\nQuant","FlashLFQ\nv1.0")


pip_plot <- ggplot(data= pip_table_subset, aes(x = PipFilter, y = pip_counts)) +
  geom_bar(stat="identity", position=position_dodge(), fill =  alpha("#0072B2", .75)) +
  facet_grid(~DonorFilter, space = "free", scales = "free") +
  ylab("Number of PIP Events") +
  xlab("PIP FDR Cutoff (%)") +
  scale_y_continuous(labels = comma) +
  theme_bw()  +
  theme(panel.grid.major.y = element_line(color = "darkgray"),
        legend.position = "none",
        axis.text = element_text(size=ax_text_size),
        axis.title=element_text(size=ax_title_size),
        plot.title = element_text(size=plot_title_size),
        strip.text.x = element_text(size = strip_text_size))

ggsave(filename = paste0(base_out_dir, "PIP_Plot_sub", ".png"),
       plot = pip_plot,
       units = "px",
       width = 2700,
       height = 2100)

# Specific Comparison FlashLFQ v2.0 - Figure 4A----
bar_height = 160
ymax = 180

get_fc_hist_specific <- function(results, comp, fc = 2.5, e_coli_count = 0, fc_delta = 0, title = "")
{
  if(length(e_coli_count) == 0)
  {
    e_coli_count <- 0
  }
  
  plot <- ggplot(data = results, aes(fill = Organism, x = Log2_fc)) +
    geom_histogram(alpha = 0.44, position = "identity", binwidth = 0.02, show.legend = T) +
    coord_cartesian(xlim = c(0.4, 1.9), ylim = c(-2, ymax)) +
    scale_fill_manual(
      name = "Species",
      values = c("Homo sapiens" = "black", "Escherichia coli (strain K12)" = "#009970"),
      labels = c("E. coli", "H. sapiens")) +
    xlab("Log2 Fold Change") +
    ylab("Number of Peptides") +
    theme_bw()  +
    theme(panel.grid.major.y = element_line(color = "darkgray"),
          legend.position = "right",
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
  #+ ggtitle(title) plot.title = element_text(hjust = 0.5, size = title_size),
  
  plot <- plot + theme(legend.text = element_text(size = 14),
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

new_flash_w_pep_path <- r"(D:\PXD003881_IonStar_SpikeIn\FlashLFQ_7772_Normed_Donor05_Pip2p5\QuantifiedPeptides.tsv)"
ymax <- 200
text_height <- 185
bar_height <- 170

conditions <- c(all_conditions[k], all_conditions[y])
pip_table <- get_peptide_table(new_flash_w_pep_path, conditions)
fv2_w_pip <- count_peps_at_fdp(pip_table, conditions, fc, order_by_fc = FALSE, sanity_check = F)
fv2_w_pip_plot <- get_fc_hist_specific(fv2_w_pip$peptide_table, comparison, fc,
                                       e_coli_count = fv2_w_pip$true_positives, fc_delta = fv2_w_pip$fc_delta,
                                       title = "FlashLFQ v2.0 with PIP - 2.5x vs 1x")

specific_out_dir <- r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\FcSpecificComparisons\)"

ggsave(filename = paste0(specific_out_dir, "fv2_w_pip", ".png"),
       plot = fv2_w_pip_plot,
       units = "px",
       width = 2100,
       height = 1500)

#Analysis of the IonStar spike-in data
library(dplyr)
library(ggplot2)
library(limma)
library(tibble)
library(matrixStats)
library(scales)

#Determines the P value threshold for the limma analysis. 0.01 For S2-A,C; 0.05 for S2-B,D
global_p_value_threshold <- 0.01

#Functions ----

calculate_DEA_results <- function(results_table, fc_threshold = 1.0, p_val_threshold = global_p_value_threshold)
{
  true_pos <- sum(results_table$adj.P.Val < p_val_threshold
                  & results_table$Organism == "Escherichia coli (strain K12)" 
                  & results_table$logFC >= log2(fc_threshold))
  
  false_pos <- sum(results_table$adj.P.Val < p_val_threshold 
                   & results_table$Organism == "Homo sapiens"
                   & results_table$logFC >= log2(fc_threshold))
  
  
  return(list(true_positive_count = true_pos,
              false_positive_count = false_pos,
              fdr = false_pos / (true_pos+false_pos)))
}

calc_DEA_2 <- function(peptide_table, fc_threshold = 1.0, fdp = global_p_value_threshold)
{
  peptide_table <- peptide_table[order(peptide_table$adj.P.Val, decreasing = F), ]
  
  peptide_table <- within(peptide_table, e_coli_cumsum <- cumsum(
    (Organism == "Escherichia coli (strain K12)") & (logFC >= log2(fc_threshold))))
  peptide_table <- within(peptide_table, human_cumsum <- cumsum(
    (Organism == "Homo sapiens"))) #& (logFC >= log2(fc_threshold))))
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
  human_count =  pruned_peptide_table$human_cumsum %>% tail(n=1)
  
  return(list(true_positive_count = e_coli_count, 
              false_positive_count = human_count))
}

LIMMA_peptide_Analysis <- function(peptides, condition_ids)
{
  groupa <- get_col_names(peptides, condition_ids[1])
  groupb <- get_col_names(peptides, condition_ids[2])
  peptide_raw <- peptides[c(groupa, groupb)]
  peptide_raw[peptide_raw == 0] <- NA
  peptide_raw <- log2(peptide_raw)
  rownames(peptide_raw) <- peptides$Sequence
  peptide_raw <- peptide_raw[rowSums(is.na(peptide_raw[groupa])) <= 2, ]
  peptide_raw <- peptide_raw[rowSums(is.na(peptide_raw[groupb])) <= 2, ]
  peptide_raw <- peptide_raw %>% as.matrix()
  design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,0,1,1,1,1))
  
  fit <- lmFit(peptide_raw, design)
  fit <- eBayes(fit, trend = T)
  results <- topTable(fit, coef=2, number = 5000)
  results$Sequence <- rownames(results)
  results <- merge(x = results, y = peptides, by = "Sequence")
  results <- results[order(results$adj.P.Val, decreasing = FALSE), ]
  row.names(results) <- NULL
  
  return(results)
}

LIMMA_protein_Analysis <- function(proteins, condition_ids = c("ecoli_D_3", "ecoli_E_3"))
{
  label_columns <- colnames(proteins)[1:3]
  
  groupa <- colnames(proteins)[grep(condition_ids[1], colnames(proteins))]
  groupa <- groupa[grep("Intensity", groupa)]
  groupb <- colnames(proteins)[grep(condition_ids[2], colnames(proteins))]
  groupb <- groupb[grep("Intensity", groupa)]
  proteins <- proteins[c(label_columns, groupa, groupb)]
  
  proteins[proteins==0] <- NA
  proteins[c(groupa, groupb)] <- log2(proteins[c(groupa, groupb)])
  protein_raw <- proteins[c(groupa, groupb)]
  rownames(protein_raw) <- proteins[,"Protein.Groups"]
  
  protein_raw <- protein_raw[rowSums(is.na(protein_raw[groupa])) <= 2, ]
  protein_raw <- protein_raw[rowSums(is.na(protein_raw[groupb])) <= 2, ]
  protein_raw <- protein_raw %>% as.matrix()
  
  design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,0,1,1,1,1))
  
  fit <- lmFit(protein_raw, design)
  fit <- eBayes(fit, trend = T)
  results <- topTable(fit, coef=2, number = 1000)
  results$`Protein.Groups` <- rownames(results)
  results <- merge(x = results, y = proteins, by = "Protein.Groups")
  results <- results[order(results$adj.P.Val, decreasing = FALSE), ]
  row.names(results) <- NULL
  
  return(results)
}

IonQuant_peptide_Analysis <- function(ions, condition_ids = c("A_", "E_"))
{
  groupa <-  colnames(ions)[grepl(condition_ids[1], colnames(ions))]
  groupb <- colnames(ions)[grepl(condition_ids[2], colnames(ions))]
  info_col <- c("Peptide.Sequence","Modified.Sequence",  "Protein",                
                "Protein.ID",  "Entry.Name", "Gene", "Protein.Description",    
                "Mapped.Genes", "Mapped.Proteins")
  
  peptide <- ions[c(info_col, groupa, groupb)]
  peptide[peptide==0] <- NA
  peptide[c(groupa, groupb)] <- log2(peptide[c(groupa, groupb)])
  peptide$Sequence <- ions$Modified.Sequence
  peptide$na_count <- rowSums(is.na(peptide))
  
  peptide <- peptide %>%
    group_by(Sequence) %>%
    slice(which.min(na_count))
  
  peptide <- peptide[rowSums(is.na(peptide[groupa])) <= 2, ]
  peptide <- peptide[rowSums(is.na(peptide[groupb])) <= 2, ]
  peptide_raw <- peptide[c(groupa, groupb)]
  rownames(peptide_raw) <- peptide$Sequence
  
  design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,0,1,1,1,1))
  
  fit <- lmFit(peptide_raw, design)
  fit <- eBayes(fit, trend = T)
  results <- topTable(fit, coef=2, number = 5000)
  results$Sequence <- rownames(results)
  results <- merge(x = results, y = peptide, by = "Sequence")
  results <- results[order(results$adj.P.Val, decreasing = FALSE), ]
  row.names(results) <- NULL
  
  return(results)
}

calculate_DEA_results_iq <- function(results_table, fc_threshold = 1.0, p_val_threshold = global_p_value_threshold)
{
  true_pos <- sum(results_table$adj.P.Val < p_val_threshold
                  & grepl("ECOLI", results_table$Protein) 
                  & results_table$logFC >= log2(fc_threshold))
  
  false_pos <- sum(results_table$adj.P.Val < p_val_threshold 
                   & grepl("HUMAN", results_table$Protein) 
                   & results_table$logFC >= log2(fc_threshold))
  
  
  return(list(true_positive_count = true_pos,
              false_positive_count = false_pos,
              fdr = false_pos / (true_pos+false_pos)))
}

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

#Define variables ----

base_path = "D:/PXD003881_IonStar_SpikeIn/FlashLFQ_7772_Normed_"
donors = c("Donor02", "Donor05","Donor1")
pips = c("_Pip1", "_Pip2p5", "_Pip5", "_Pip100")
ion_conditions <- c("A_", "B_", "C_",
                    "D_", "E_")

base_out_dir <- r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\Limma_Analysis_10_25\)"
dir.create(base_out_dir, showWarnings = FALSE)

# IonQuant Protein Analyis ----

IonQuant_protein_Analysis <- function(iq_protein_path, condition_ids = c("A_", "E_"))
{
  iq_pro <- read.csv(iq_protein_path, sep = '\t')
  iq_pro <- subset(iq_pro, Organism %in% c("Escherichia coli (strain K12)","Homo sapiens"))
  
  info_col_pro <- c("Protein",                
                    "Protein.ID",  "Entry.Name", "Gene", "Organism" )
  intensity_col <- colnames(iq_pro)[grepl("Intensity",colnames(iq_pro))]
  intensity_col <- intensity_col[grepl("MaxLFQ", intensity_col)]
  ions <- iq_pro[c(info_col_pro, intensity_col)]
  
  groupa <-  colnames(ions)[grepl(condition_ids[1], colnames(ions))]
  groupb <- colnames(ions)[grepl(condition_ids[2], colnames(ions))]
  info_col_pro <- c("Protein",                
                    "Protein.ID",  "Entry.Name", "Gene" )
  
  peptide <- ions[c(info_col_pro, groupa, groupb)]
  peptide[peptide==0] <- NA
  peptide[c(groupa, groupb)] <- log2(peptide[c(groupa, groupb)])
  peptide$na_count <- rowSums(is.na(peptide))
  
  peptide <- peptide[rowSums(is.na(peptide[groupa])) <= 2, ]
  peptide <- peptide[rowSums(is.na(peptide[groupb])) <= 2, ]
  peptide_raw <- peptide[c(groupa, groupb)]
  rownames(peptide_raw) <- peptide$Protein %>% as.matrix()
  
  design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,0,1,1,1,1))
  
  fit <- lmFit(peptide_raw, design)
  fit <- eBayes(fit, trend = T)
  results <- topTable(fit, coef=2, number = 2000)
  results$Protein <- rownames(results)
  results <- merge(x = results, y = peptide, by = "Protein")
  results <- results[order(results$adj.P.Val, decreasing = FALSE), ]
  row.names(results) <- NULL
  
  return(results)
}

base_path_iq = "D:/PXD003881_IonStar_SpikeIn/IonQuant_Normed"
pips = c("_NoPip", "_Pip1", "_Pip2p5", "_Pip5", "_Pip100")
file_path_iq = "/combined_protein.tsv"

tp_counts_iq <- c()
fp_counts_iq <- c()
fdp_iq <- c()
pip_settings_iq <- c()

for(j in 1:5)
{
  ion_path = paste0(base_path_iq, pips[j], file_path_iq)
  
  true_pos_running_total = 0
  false_pos_running_total = 0
  for(k in 1:4)
  {
    for(y in (k+1):5)
    {
      iq_pro_res <- IonQuant_protein_Analysis(ion_path, c(ion_conditions[k], ion_conditions[y]))
      dea_counts = calculate_DEA_results_iq(iq_pro_res)
      true_pos_running_total = true_pos_running_total + dea_counts$true_positive_count
      false_pos_running_total = false_pos_running_total + dea_counts$false_positive_count
    }
  }
  
  fdr_iq_pro <- false_pos_running_total / (true_pos_running_total + false_pos_running_total)
  
  pip_settings_iq <- c(pip_settings_iq, pips[j])
  tp_counts_iq <- c(tp_counts_iq, true_pos_running_total)
  fp_counts_iq <- c(fp_counts_iq, false_pos_running_total)
  fdp_iq <- c(fdp_iq, fdr_iq_pro)
  print(fdr_iq_pro)
}

# FlashLFQ + PIP-ECHO Protein Analysis ----

all_conditions <- c("ecoli_A_3", "ecoli_B_3", "ecoli_C_3",
                    "ecoli_D_3", "ecoli_E_3")

base_path = "D:/PXD003881_IonStar_SpikeIn/FlashLFQ_7772_Normed_"
donors = c("Donor02", "Donor05","Donor1")
pips = c("_Pip1", "_Pip2p5", "_Pip5", "_Pip100")

protein_file_path = "/QuantifiedProteins.tsv"

pip_settings = c()
true_positives_pro = c()
false_positives_pro = c()
fdr_pro = c()
for(i in 1:3)
{
  for(j in 1:4)
  {
    protein_path =paste0(base_path, donors[i], pips[j], protein_file_path)
    pip_settings <- c(pip_settings, paste0(donors[i], pips[j]))
    true_pos_running_total = 0
    false_pos_running_total = 0
    
    proteins <- read.table(protein_path, header = TRUE, sep = '\t', fill = TRUE)
    if("Protein.Accession" %in% colnames(proteins))
    {
      names(proteins)[names(proteins) == "Protein.Accession"] <- "Protein.Groups"
    }
    proteins <- proteins[!grepl("DECOY", proteins$Protein.Groups),]
    proteins <- subset(proteins, Organism %in%
                         c("Escherichia coli (strain K12)",
                           "Homo sapiens"))
    for(k in 1:4)
    {
      for(y in (k+1):5)
      {
        limma_protein_results <- LIMMA_protein_Analysis(proteins,
                                                        c(all_conditions[k], all_conditions[y]))
        dea_counts = calculate_DEA_results(limma_protein_results)
        new_dea_counts = calc_DEA_2(limma_protein_results)
        true_pos_running_total = true_pos_running_total + dea_counts$true_positive_count
        false_pos_running_total = false_pos_running_total + dea_counts$false_positive_count
      }
    }
    
    true_positives_pro <- c(true_positives_pro, true_pos_running_total)
    false_positives_pro <- c(false_positives_pro, false_pos_running_total)
    fdr_pro <- c(fdr_pro, false_pos_running_total / (false_pos_running_total + true_pos_running_total))
    
    print(paste0(paste0(donors[i], pips[j]), " - True Positives = ", true_pos_running_total))
  }
}

# FlashLFQ v1.0 Protein Analysis ----

flash_v1_base <- r"(D:\PXD003881_IonStar_SpikeIn\FlashLFQ_Current_Normed_NewPepQ)"
pips = c("_NoMBR", "")
file_path = "/QuantifiedProteins.tsv"
pip_settings_v1 = c()
tp_counts_v1 = c()
fp_counts_v1 = c()
fdp_v1 = c()

for(i in 1:2)
{
  true_pos_running_total = 0
  false_pos_running_total = 0
  pip_settings_v1 <- c(pip_settings_v1, pips[i])
  old_flash_path = paste0(flash_v1_base, pips[i], file_path)
  proteins <- read.table(old_flash_path, header = TRUE, sep = '\t', fill = TRUE)
  if("Protein.Accession" %in% colnames(proteins))
  {
    names(proteins)[names(proteins) == "Protein.Accession"] <- "Protein.Groups" 
  }
  proteins <- proteins[!grepl("DECOY", proteins$Protein.Groups),]
  proteins <- subset(proteins, Organism %in% 
                       c("Escherichia coli (strain K12)", 
                         "Homo sapiens"))
  for(k in 1:4)
  {
    for(y in (k+1):5)
    {
      limma_protein_results <- LIMMA_protein_Analysis(proteins,
                                                      c(all_conditions[k], all_conditions[y]))
      dea_counts = calculate_DEA_results(limma_protein_results)
      true_pos_running_total = true_pos_running_total + dea_counts$true_positive_count
      false_pos_running_total = false_pos_running_total + dea_counts$false_positive_count
    }
  }
  tp_counts_v1 <- c(tp_counts_v1, true_pos_running_total)
  fp_counts_v1 <- c(fp_counts_v1, false_pos_running_total)
  fdp_v1 <- c(fdp_v1, (false_pos_running_total / (false_pos_running_total + true_pos_running_total)))
}

# Combine Protein Results ----

results_list <- list(pip_condition = c(paste0("FlashLFQ v1.0", pip_settings_v1),
                                       paste0("IonQuant", pip_settings_iq),
                                       pip_settings),
                     true_positives = c(tp_counts_v1, tp_counts_iq, true_positives_pro),
                     false_positives = c(fp_counts_v1, fp_counts_iq, false_positives_pro),
                     fdr = c(fdp_v1, fdp_iq, fdr_pro))
combined_results <- as.data.frame(results_list)
combined_results$pip_condition <- factor(combined_results$pip_condition,
                                         levels = c(paste0("FlashLFQ v1.0", pip_settings_v1),
                                                    paste0("IonQuant", pip_settings_iq),
                                                    pip_settings))
combined_results$DonorFilter <- factor(c(rep("FlashLFQ\nv1.0", 2), rep("IonQuant", 5),
                                         rep("Donor FDR Cutoff = 0.2%", 4),
                                         rep("Donor FDR Cutoff = 0.5%", 4),
                                         rep("Donor FDR Cutoff = 1%", 4)),
                                       levels = c("FlashLFQ\nv1.0", "IonQuant",
                                                  "Donor FDR Cutoff = 0.2%",
                                                  "Donor FDR Cutoff = 0.5%",
                                                  "Donor FDR Cutoff = 1%"))

combined_results$PipFilter <- factor(c("No PIP", 100, "No PIP", 1, 2.5, 5, 100, 1, 2.5, 5, 100, 1, 2.5, 5, 100, 1, 2.5, 5, 100),
                                     levels = c("No PIP", 1, 2.5, 5, 100))


# Select the condition we want to include in the main figure
mask <- (combined_results$PipFilter == 1 & combined_results$DonorFilter == "Donor FDR Cutoff = 0.2%") |
  (combined_results$PipFilter == 2.5 & combined_results$DonorFilter == "Donor FDR Cutoff = 0.5%") |
  (combined_results$PipFilter == 5 & combined_results$DonorFilter == "Donor FDR Cutoff = 1%") |
  (combined_results$PipFilter == 100 & combined_results$DonorFilter == "Donor FDR Cutoff = 1%") |
  combined_results$DonorFilter == "FlashLFQ\nv1.0" | combined_results$DonorFilter == "IonQuant" 

combined_results_subset <- combined_results[mask,]

# rename and re-order 
#levels(combined_results_subset$DonorFilter)[levels(combined_results_subset$DonorFilter) == "FlashLFQ v1.0"] <- "FlashLFQ\nv1.0"
levels(combined_results_subset$DonorFilter)[levels(combined_results_subset$DonorFilter) %in% 
                                              c("Donor FDR Cutoff = 0.2%", 
                                                "Donor FDR Cutoff = 0.5%",
                                                "Donor FDR Cutoff = 1%")] <- "FlashLFQ +\nPIP-ECHO"

combined_results_subset$DonorFilter <- factor(combined_results_subset$DonorFilter,
                                              levels = c("FlashLFQ +\nPIP-ECHO",
                                                         "IonQuant",
                                                         "FlashLFQ\nv1.0"))
flash_no_pip <- combined_results_subset[combined_results_subset$DonorFilter == "FlashLFQ\nv1.0" 
                                        & combined_results_subset$PipFilter == "No PIP",]
flash_no_pip$DonorFilter <- "FlashLFQ +\nPIP-ECHO"
combined_results_subset <- rbind(combined_results_subset, flash_no_pip)

# IonQuant Peptide Analysis ----

IonQuant_peptide_Analysis <- function(ion_path, condition_ids = c("A_", "E_"))
{
  ions <- read.csv(ion_path, sep = '\t')
  info_col <- c("Peptide.Sequence","Modified.Sequence",  "Protein",                
                "Protein.ID",  "Entry.Name", "Gene", "Protein.Description",    
                "Mapped.Genes", "Mapped.Proteins"  )
  intensity_col <- colnames(ions)[grepl("Intensity",colnames(ions))]
  ions <- ions[c(info_col, intensity_col)]
  ions <- ions[grepl("HUMAN",ions$Entry.Name) | grepl("ECOLI", ions$Entry.Name),]
  
  groupa <-  colnames(ions)[grepl(condition_ids[1], colnames(ions))]
  groupb <- colnames(ions)[grepl(condition_ids[2], colnames(ions))]
  info_col <- c("Peptide.Sequence","Modified.Sequence",  "Protein",                
                "Protein.ID",  "Entry.Name", "Gene", "Protein.Description",    
                "Mapped.Genes", "Mapped.Proteins")
  
  peptide <- ions[c(info_col, groupa, groupb)]
  peptide[peptide==0] <- NA
  peptide[c(groupa, groupb)] <- log2(peptide[c(groupa, groupb)])
  peptide$Sequence <- ions$Modified.Sequence
  peptide$na_count <- rowSums(is.na(peptide))
  
  peptide <- peptide %>%
    group_by(Sequence) %>%
    slice(which.min(na_count))
  
  peptide <- peptide[rowSums(is.na(peptide[groupa])) <= 2, ]
  peptide <- peptide[rowSums(is.na(peptide[groupb])) <= 2, ]
  peptide_raw <-peptide[c(groupa, groupb)] %>% as.matrix()
  rownames(peptide_raw) <- peptide$Sequence
  
  design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,0,1,1,1,1))
  
  fit <- lmFit(peptide_raw, design)
  fit <- eBayes(fit, trend = T)
  results <- topTable(fit, coef=2, number = 5000)
  results$Sequence <- rownames(results)
  results <- merge(x = results, y = peptide, by = "Sequence")
  results <- results[order(results$adj.P.Val, decreasing = FALSE), ]
  
  if(!is.null(results) && nrow(results) > 0) {
    row.names(results) <- NULL}
  
  return(list("results" = results, "quant_count" = nrow(peptide)))
}

base_path_iq = "D:/PXD003881_IonStar_SpikeIn/IonQuant_Normed"
pips = c("_NoPip", "_Pip1", "_Pip2p5", "_Pip5", "_Pip100")
ion_conditions <- c("A_", "B_", "C_",
                    "D_", "E_")
file_path_iq = "/combined_ion.tsv"

tp_counts_iq_pep <- c()
fp_counts_iq_pep <- c()
fdp_iq_pep <- c()
quant_count_iq_pep <- c()
pip_settings_iq_pep <- c()

for(j in 1:5)
{
  ion_path = paste0(base_path_iq, pips[j], file_path_iq)
  
  true_pos_running_total = 0
  false_pos_running_total = 0
  total_peptides_quantified_rt = 0
  for(k in 1:4)
  {
    for(y in (k+1):5)
    {
      iq_pro_res <- IonQuant_peptide_Analysis(ion_path, c(ion_conditions[k], ion_conditions[y]))
      dea_counts = calculate_DEA_results_iq(iq_pro_res$results)
      true_pos_running_total = true_pos_running_total + dea_counts$true_positive_count
      false_pos_running_total = false_pos_running_total + dea_counts$false_positive_count
      total_peptides_quantified_rt = total_peptides_quantified_rt + iq_pro_res$quant_count
    }
  }
  
  fdr_iq_pro <- false_pos_running_total / (true_pos_running_total + false_pos_running_total)
  
  quant_count_iq_pep <- c(quant_count_iq_pep, total_peptides_quantified_rt)
  pip_settings_iq_pep <- c(pip_settings_iq_pep, pips[j])
  tp_counts_iq_pep <- c(tp_counts_iq_pep, true_pos_running_total)
  fp_counts_iq_pep <- c(fp_counts_iq_pep, false_pos_running_total)
  fdp_iq_pep <- c(fdp_iq_pep, fdr_iq_pro)
  print(fdr_iq_pro)
}

# FlashLFQ + PIP-Echo Peptide Analysis ----

LIMMA_peptide_Analysis <- function(peptides, condition_ids)
{
  groupa <- get_col_names(peptides, condition_ids[1])
  groupb <- get_col_names(peptides, condition_ids[2])
  peptide_raw <- peptides[c(groupa, groupb)]
  peptide_raw[peptide_raw == 0] <- NA
  peptide_raw <- log2(peptide_raw)
  rownames(peptide_raw) <- peptides$Sequence
  peptide_raw <- peptide_raw[rowSums(is.na(peptide_raw[groupa])) <= 2, ]
  peptide_raw <- peptide_raw[rowSums(is.na(peptide_raw[groupb])) <= 2, ]
  peptide_raw <- peptide_raw %>% as.matrix()
  design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,0,1,1,1,1))
  
  fit <- lmFit(peptide_raw, design)
  fit <- eBayes(fit, trend = T)
  results <- topTable(fit, coef=2, number = 5000)
  results$Sequence <- rownames(results)
  results <- merge(x = results, y = peptides, by = "Sequence")
  results <- results[order(results$adj.P.Val, decreasing = FALSE), ]
  row.names(results) <- NULL
  
  return(list("results" = results, "quant_count" = nrow(peptide_raw)))
}

all_conditions <- c("ecoli_A_3", "ecoli_B_3", "ecoli_C_3",
                    "ecoli_D_3", "ecoli_E_3")

base_path = "D:/PXD003881_IonStar_SpikeIn/FlashLFQ_7772_Normed_"
donors = c("Donor02", "Donor05","Donor1")
pips = c("_Pip1", "_Pip2p5", "_Pip5", "_Pip100")
peptide_file = "/QuantifiedPeptides.tsv"

pip_settings = c()
true_positives = c()
false_positives = c()
quant_count = c()
fdr_pep = c()

for(i in 1:3)
{
  for(j in 1:4)
  {
    peptide_path = paste0(base_path, donors[i], pips[j], peptide_file)
    
    peptides <- read.table(peptide_path, header = TRUE, sep = '\t')
    peptides <- peptides[!grepl("DECOY", peptides$Protein.Groups),]
    peptides <- subset(peptides, Organism %in% 
                         c("Escherichia coli (strain K12)", 
                           "Homo sapiens"))
    label_columns <- colnames(peptides)[c(1,3,4,5)]
    intensity_columns <- colnames(peptides)[grepl("Intensity", colnames(peptides))]
    peptides <- peptides[c(label_columns, intensity_columns)]
    
    pip_settings <- c(pip_settings, paste0(donors[i], pips[j]))
    true_pos_running_total = 0
    false_pos_running_total = 0
    total_peptides_quantified_rt = 0
    for(k in 1:4)
    {
      for(y in (k+1):5)
      {
        limma_pep_results <- LIMMA_peptide_Analysis(peptides,
                                                    c(all_conditions[k], all_conditions[y]))
        dea_counts = calculate_DEA_results(limma_pep_results$results)
        new_dea_counts = calc_DEA_2(limma_pep_results$results)
        true_pos_running_total = true_pos_running_total + dea_counts$true_positive_count
        false_pos_running_total = false_pos_running_total + dea_counts$false_positive_count
        
        total_peptides_quantified_rt = total_peptides_quantified_rt + limma_pep_results$quant_count
      }
    }
    
    true_positives <- c(true_positives, true_pos_running_total)
    false_positives <- c(false_positives, false_pos_running_total)
    fdr_pep <- c(fdr_pep, false_pos_running_total / (false_pos_running_total + true_pos_running_total))
    quant_count <- c(quant_count, total_peptides_quantified_rt)
    
    print(paste0(paste0(donors[i], pips[j]), " - True Positives = ", true_pos_running_total))
  }
}


# FlashLFQ v1.0 Peptide Analysis ----

flash_v1_base <- r"(D:\PXD003881_IonStar_SpikeIn\FlashLFQ_Current_Normed_NewPepQ)"
pips = c("_NoMBR", "")
file_path = "/QuantifiedPeptides.tsv"
pip_settings_v1 = c()
tp_counts_v1 = c()
fp_counts_v1 = c()
fdp_v1 = c()

for(i in 1:2)
{
  true_pos_running_total = 0
  false_pos_running_total = 0
  pip_settings_v1 <- c(pip_settings_v1, pips[i])
  old_flash_path = paste0(flash_v1_base, pips[i], file_path)
  peptides <- read.table(old_flash_path, header = TRUE, sep = '\t', fill = TRUE)
  if("Protein.Accession" %in% colnames(peptides))
  {
    names(peptides)[names(peptides) == "Protein.Accession"] <- "Protein.Groups" 
  }
  peptides <- peptides[!grepl("DECOY", peptides$Protein.Groups),]
  peptides <- subset(peptides, Organism %in% 
                       c("Escherichia coli (strain K12)", 
                         "Homo sapiens"))
  label_columns <- colnames(peptides)[c(1,3,4,5)]
  intensity_columns <- colnames(peptides)[grepl("Intensity", colnames(peptides))]
  peptides <- peptides[c(label_columns, intensity_columns)]
  for(k in 1:4)
  {
    for(y in (k+1):5)
    {
      limma_pep_results <- LIMMA_peptide_Analysis(peptides,
                                                  c(all_conditions[k], all_conditions[y]))
      dea_counts = calculate_DEA_results(limma_pep_results$results)
      true_pos_running_total = true_pos_running_total + dea_counts$true_positive_count
      false_pos_running_total = false_pos_running_total + dea_counts$false_positive_count
    }
  }
  
  tp_counts_v1 <- c(tp_counts_v1, true_pos_running_total)
  fp_counts_v1 <- c(fp_counts_v1, false_pos_running_total)
  fdp_v1 <- c(fdp_v1, (false_pos_running_total / (false_pos_running_total + true_pos_running_total)))
}

# Combine Peptide Results ----
pep_results_list <- list(pip_condition = c(paste0("FlashLFQ v1.0", pip_settings_v1),
                                           paste0("IonQuant", pip_settings_iq_pep),
                                           pip_settings),
                         true_positives = c(tp_counts_v1, tp_counts_iq_pep, true_positives),
                         false_positives = c(fp_counts_v1, fp_counts_iq_pep, false_positives),
                         fdr = c(fdp_v1, fdp_iq_pep, fdr_pep))
pep_results <- as.data.frame(pep_results_list)
pep_results$pip_condition <- factor(pep_results$pip_condition,
                                    levels = c(paste0("FlashLFQ\nv1.0", pip_settings_v1),
                                               paste0("IonQuant", pip_settings_iq_pep),
                                               pip_settings))
pep_results$DonorFilter <- factor(c(rep("FlashLFQ\nv1.0", 2), rep("IonQuant", 5),
                                    rep("Donor FDR Cutoff = 0.2%", 4),
                                    rep("Donor FDR Cutoff = 0.5%", 4),
                                    rep("Donor FDR Cutoff = 1%", 4)),
                                  levels = c("FlashLFQ\nv1.0", "IonQuant",
                                             "Donor FDR Cutoff = 0.2%",
                                             "Donor FDR Cutoff = 0.5%",
                                             "Donor FDR Cutoff = 1%"))

pep_results$PipFilter <- factor(c("No PIP", 100, "No PIP", 1, 2.5, 5, 100, 1, 2.5, 5, 100, 1, 2.5, 5, 100, 1, 2.5, 5, 100),
                                levels = c("No PIP", 1, 2.5, 5, 100))

# Select the condition we want to include in the main figure
mask <- (pep_results$PipFilter == 1 & pep_results$DonorFilter == "Donor FDR Cutoff = 0.2%") |
  (pep_results$PipFilter == 2.5 & pep_results$DonorFilter == "Donor FDR Cutoff = 0.5%") |
  (pep_results$PipFilter == 5 & pep_results$DonorFilter == "Donor FDR Cutoff = 1%") |
  (pep_results$PipFilter == 100 & pep_results$DonorFilter == "Donor FDR Cutoff = 1%") |
  pep_results$DonorFilter == "FlashLFQ\nv1.0" | pep_results$DonorFilter == "IonQuant" 

pep_results_subset <- pep_results[mask,]

# rename and re-order 
levels(pep_results_subset$DonorFilter)[levels(pep_results_subset$DonorFilter) %in% 
                                         c("Donor FDR Cutoff = 0.2%", 
                                           "Donor FDR Cutoff = 0.5%",
                                           "Donor FDR Cutoff = 1%")] <- "FlashLFQ +\nPIP-ECHO"

pep_results_subset$DonorFilter <- factor(pep_results_subset$DonorFilter,
                                         levels = c("FlashLFQ +\nPIP-ECHO",
                                                    "IonQuant",
                                                    "FlashLFQ\nv1.0"))
flash_no_pip <- pep_results_subset[pep_results_subset$DonorFilter == "FlashLFQ\nv1.0" 
                                   & pep_results_subset$PipFilter == "No PIP",]
flash_no_pip$DonorFilter <- "FlashLFQ +\nPIP-ECHO"
pep_results_subset <- rbind(pep_results_subset, flash_no_pip)

# All plots ----

legend_text_size = 10
ax_text_size = 12
strip_text_size = 13
ax_title_size = 16
plot_title_size = 18

fdp_text_size = 4

# Proteins
crs_stash <- combined_results_subset
crs_stash$PipFilter <- factor(c("No PIP", "N/A", "No PIP", 1, 2.5, 5, "N/A", 1, 2.5, 5, "N/A", "No PIP"),
                              levels = c("No PIP", 1, 2.5, 5, "N/A"))


base_out_dir <- r"(C:\Users\Alex\Source\Repos\MBR_metamorpheus\figures\MockFigures\R_Plots\Limma_Analysis_10_27\)"
dir.create(base_out_dir, showWarnings = FALSE)


crs_stash <- combined_results_subset
crs_stash$PipFilter <- factor(c("No PIP", "N/A", "No PIP", 1, 2.5, 5, "N/A", 1, 2.5, 5, "N/A", "No PIP"),
                              levels = c("No PIP", 1, 2.5, 5, "N/A"))

crs_plot <- ggplot(data= crs_stash, aes(x = PipFilter, y = true_positives, fill = 100*fdr)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(x = PipFilter, y = true_positives, label = paste0(format(fdr*100, digits = 1), "%")),
            position=position_dodge(width = 1),
            vjust = -0.5,
            size = fdp_text_size) +
  facet_grid(~DonorFilter, space = "free", scales = "free") +
  scale_fill_gradient(high = "red", low= "gray", name="False\nDiscovery\nProportion (%)") +
  ylab("Number of differentially abundant E. coli proteins" ) +
  xlab("PIP FDR Cutoff (%)") +
  ggtitle("Differentially Abundant Proteins, FDR = 0.05") +
  coord_cartesian(ylim=c(2200, 6400)) +
  #coord_cartesian(ylim=c(2200, 5500)) +
  scale_y_continuous(labels = comma) +
  theme_bw()   +
  theme(axis.text = element_text(size=ax_text_size),
        axis.title=element_text(size=ax_title_size),
        plot.title = element_text(size=plot_title_size),
        strip.text.x = element_text(size = strip_text_size),
        legend.title=element_text(size=ax_text_size),
        legend.text=element_text(size= legend_text_size))

ggsave(filename = paste0(base_out_dir, "Limma_Protein_Level_05", ".png"),
       plot = crs_plot,
       units = "px",
       width = 2400,
       height = 1900)


# Peptides

prs_stash <- pep_results_subset
prs_stash$PipFilter <- factor(c("No PIP", "N/A", "No PIP", 1, 2.5, 5, "N/A", 1, 2.5, 5, "N/A", "No PIP"),
                              levels = c("No PIP", 1, 2.5, 5, "N/A"))

prs_plot <- ggplot(data= prs_stash, aes(x = PipFilter, y = true_positives, fill = fdr*100)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(x = PipFilter, y = true_positives, label = paste0(format(fdr*100, digits = 2), "%")),
            position=position_dodge(width = 1),
            vjust = -0.5,
            size = fdp_text_size) +
  facet_grid(~DonorFilter, space = "free", scales = "free") +
  scale_fill_gradient(high = "red", low= "gray", name="False\nDiscovery\nProportion (%)") +
  ylab("Number of differentially abundant E. coli peptides" ) +
  xlab("PIP FDR Cutoff (%)") +
  ggtitle("Differentially Abundant Peptides, FDR = 0.05") +
  coord_cartesian(ylim=c(6000, 19000)) +
  scale_y_continuous(labels = comma) +
  theme_bw()   +
  theme(axis.text = element_text(size=ax_text_size),
        axis.title=element_text(size=ax_title_size),
        plot.title = element_text(size=plot_title_size),
        strip.text.x = element_text(size = strip_text_size),
        legend.title=element_text(size=ax_text_size),
        legend.text=element_text(size= legend_text_size))


ggsave(filename = paste0(base_out_dir, "Limma_Peptide_Level_05", ".png"),
       plot = prs_plot,
       units = "px",
       width = 2400,
       height = 1900)
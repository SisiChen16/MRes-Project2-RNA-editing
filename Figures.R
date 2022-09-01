# Load packages
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(dplyr)
library(fst)
library(broom)
library(ggpubr)
library(lme4)
library(stringr)
library(lmerTest)
library(pheatmap)
library(ggrepel)
library(survival)
library(survminer)

# Set theme and colors
my_theme <- theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14),
                  axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14),
                  plot.title = element_text(size = 16), legend.text = element_text(size = 12), 
                  legend.title = element_text(size =12)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_blank())

my_theme2 <- theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14),
                   axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14),
                   plot.title = element_text(size = 16), legend.text = element_text(size = 12), 
                   legend.title = element_text(size =12)) 

hot_cold_color = c("#F8766D", "#00BFC4")
normal_tumour_color = c("#006837","#FF7F00")
pri_meta_color = c("#993404","#1F78B4")
lohhla_color = c("#A6CEE3","#1F78B4")
variant_present_color = c("#B15928","#FFFF99")


# Import data 
transcript <- read.fst("2022-04-28_rsem_tpm_mat.fst")
neoantigen <- read.fst("2022-05-25_neoantigens_editing.fst")
load("merged_data.RData")
load("lohhla_data.RData")
load("neoantigen2.RData")
load("region_purity.RData")
load("neo_transcript3.RData")

danaher_scores <- read.fst("2022-06-30_danaher_scores.fst")
danaher_scores <- danaher_scores %>% select(cell_type_title, score, residual, sample_name)
colnames(danaher_scores)[4] <- "patient_region"

#---------------------------------------------------------------------------------------------------------------------
## Figure 1 General overview of RNA editing 
  # Compare RNA editing pattern by histology -------------------------------------------------
# Total RNA editing variants per region
merged_data %>% filter(!is.na(tumour_id), !is.na(histology)) %>% group_by(patient_region, histology) %>% 
  summarize(total_rna = n_distinct(rna_ids)) %>% ggplot(aes(x = histology, y = total_rna, fill = histology)) + 
  geom_boxplot(color = "grey40") + stat_compare_means() + theme_bw() + scale_fill_manual(values = Histology_cols) +  
  labs(x= "Histology", y = "Total RNA editing variants") + my_theme2

# Mean editing level 
merged_data %>% filter(all_cov > 60, effect == "nonsynonymous",!is.na(tumour_id), !is.na(histology)) %>% group_by(patient_region, histology) %>% 
  summarize(mean_editing_level = mean(all_vaf)) %>% ggplot(aes(x = histology, y = mean_editing_level, fill = histology)) + 
  geom_boxplot(color = "grey40") + stat_compare_means() + scale_fill_manual(values = Histology_cols) +  
  labs(x= "Histology", y = "Mean editing level") + theme_bw() + my_theme2
 
# Mean editing level per tumour by tumour-specific vs shared
merged_data %>% filter(all_cov > 60, !is.na(tumour_id), !is.na(histology)) %>% group_by(patient_region, histology, feature2) %>% 
  summarize(mean_editing_level = mean(all_vaf)) %>% ggplot(aes(x = feature2, y = mean_editing_level, fill = histology)) + 
  geom_boxplot(color = "grey40") + stat_compare_means() + scale_fill_manual(values = Histology_cols) +  
  labs(x= "Histology", y = "Mean editing level") + theme_bw()

# Total RNA editing variants - tumour-specific vs shared
merged_data %>% filter(effect == "nonsynonymous", !is.na(tumour_id), !is.na(histology)) %>% group_by(patient_region, histology, feature2) %>% 
  summarize(total_rna = n_distinct(rna_ids)) %>% ggplot(aes(x = feature2, y = total_rna, fill = histology )) + 
  geom_boxplot(color = "grey40") + stat_compare_means() + theme_bw() + scale_fill_manual(values = Histology_cols) + my_theme2 + 
  labs(x= "", y = "Nonsynonymous RNA editing variants")

# Correlate ADAR and APOBEC expression with RNA editing variants ---------------------------------------
ADAR_APO_exp <- transcript %>% filter(gene_id == "ADAR"| gene_id == "APOBEC3A") %>% column_to_rownames(var = "gene_id"  )
ADAR_APO_exp <- as.data.frame(t(ADAR_APO_exp))

AG_variants <- merged_data %>% filter(muttype_nostrand == "A>G") %>% group_by(patient_region) %>% summarize(nb = n_distinct(rna_ids)) %>% column_to_rownames(var = "patient_region")
AG_proportion <- merged_data %>% group_by(patient_region) %>% mutate(total_rna = n_distinct(rna_ids)) %>% filter(muttype_nostrand == "A>G") %>% group_by(patient_region) %>% summarize(AG_proportion = n_distinct(rna_ids)/total_rna) %>% distinct() %>%  column_to_rownames(var = "patient_region")
CT_proportion <- merged_data %>% group_by(patient_region) %>% mutate(total_rna = n_distinct(rna_ids)) %>% filter(muttype_nostrand == "C>T") %>% group_by(patient_region) %>% summarize(CT_proportion = n_distinct(rna_ids)/total_rna) %>% distinct() %>% column_to_rownames(var = "patient_region") 
total_rna <- merged_data %>% filter(!is.na(muttype_nostrand)) %>% group_by(patient_region) %>% summarize(total_rna = n_distinct(rna_ids)) %>% column_to_rownames(var = "patient_region")

AG_merged <- merge(AG_variants, ADAR_APO_exp, by = 0)
CT_merged <- merge(CT_proportion, ADAR_APO_exp, by = 0)
total_merged <- merge(total_rna, ADAR_APO_exp, by = 0)

# Plot the linear regression 
AG_merged %>% ggplot(aes(x = nb, y = ADAR)) + geom_point() + geom_smooth(method = "lm") + stat_cor(method = "spearman") + my_theme + labs(x = "No. of A>G RNA editing variants per region", y = "ADAR expression")
CT_merged %>% ggplot(aes(x = CT_proportion, y = APOBEC3A)) + geom_point() + geom_smooth(method = "lm") + stat_cor(method = "spearman") + my_theme + labs(x = "Proportion of C>T RNA editing variants per region", y = "APOBEC3A expression")
total_merged %>% ggplot(aes(x = total_rna, y = ADAR)) + geom_point() + geom_smooth(method = "lm") + stat_cor(method = "spearman") + my_theme + labs(x = "Total RNA editing variants per region", y = "ADAR expression")

# ** Previous exploration -- APOBEC3A expression stands out 
ADAR_APO_exp <- transcript[grep("APOBEC|ADAR", transcript$gene_id),] 
rownames(ADAR_APO_exp) <- ADAR_APO_exp$gene_id
ADAR_APO_exp$gene_id <- NULL
ADAR_APO_exp <- as.data.frame(t(ADAR_APO_exp))

AG_merged <- merge(AG_variants, ADAR_APO_exp, by = 0)
CT_merged <- merge(CT_proportion, ADAR_APO_exp, by = 0)
total_merged <- merge(total_rna, ADAR_APO_exp, by = 0)

total_merged <- total_merged %>% pivot_longer(3:17, names_to = "gene_id", values_to = "expression")
total_merged %>% ggplot(aes(x = total_rna, y = expression)) + geom_point() + geom_smooth(method = "lm") + stat_cor(method = "spearman") + my_theme2 + theme_bw() + labs(x = "Total RNA editing variants per region", y = "Gene expression") + facet_wrap(~ gene_id)

AG_merged <- AG_merged %>% pivot_longer(3:17, names_to = "gene_id", values_to = "expression")
AG_merged %>% ggplot(aes(x = nb, y = expression)) + geom_point() + geom_smooth(method = "lm") + stat_cor(method = "spearman") + my_theme2 + theme_bw() + labs(x = "No. of A>G RNA editing variants per region", y = "Gene expression")+ facet_wrap(~ gene_id)

CT_merged <- CT_merged %>% pivot_longer(3:17, names_to = "gene_id", values_to = "expression")
CT_merged %>% ggplot(aes(x = CT_proportion, y = expression)) + geom_point() + geom_smooth(method = "lm") + stat_cor(method = "spearman") + my_theme2 + theme_bw() + labs(x = "Proportion of C>T RNA editing variants per region", y = "Gene expression") + facet_wrap(~ gene_id)

# Figure 2 Identification of neo-antigenic RNA editing variants ----------------------------------------
## Pie chart -- showing the effect of editing 

data <- merged_data %>% select(rna_ids, effect) %>% filter(!is.na(effect)) %>% distinct() %>% 
  group_by(effect) %>% summarize(no_rna_editing = n_distinct(rna_ids)) %>% mutate(total = sum(no_rna_editing)) %>% mutate(perc = no_rna_editing/total) %>% arrange(perc) %>% mutate(labels = scales::percent(perc)) 

data2 <- data %>% mutate(csum = rev(cumsum(rev(perc*100))), 
                         pos = perc*100/2 + lead(csum, 1),
                         pos = if_else(is.na(pos), perc*100/2, pos))

ggplot(data, aes(x = "" , y = perc*100, fill = fct_inorder(effect))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Pastel1") +
  geom_label_repel(data = data2,
                   aes(y = pos, label = labels ),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Effect of RNA editing")) +
  theme_void() + my_theme2

## Pie chart  # filter by non-synonymous
data <- merged_data %>% select(rna_ids, muttype_nostrand) %>% filter(!is.na(muttype_nostrand)) %>% distinct() %>% group_by(muttype_nostrand) %>% summarize(no_rna_editing = n_distinct(rna_ids)) %>% mutate(total = sum(no_rna_editing)) %>% mutate(perc = no_rna_editing/total) %>% arrange(perc) %>% mutate(labels = scales::percent(perc)) 

data2 <- data %>% mutate(csum = rev(cumsum(rev(perc*100))), 
                         pos = perc*100/2 + lead(csum, 1),
                         pos = if_else(is.na(pos), perc*100/2, pos))

ggplot(data, aes(x = "" , y = perc*100, fill = fct_inorder(muttype_nostrand))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  geom_label_repel(data = data2,
                   aes(y = pos, label = labels ),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "RNA editing variant type")) +
  theme_void() + scale_fill_brewer(palette="Set3")

## Neo-antigenic RNA editing variants -- Barplot ----------------------------------------------------------------
# Barplot 
test <- merged_data %>% filter(neoantigenic == TRUE) %>% select(rna_ids, hla_allele_type, tumour_id, all_vaf, EL_score, symbol,combine) %>% distinct() %>% group_by(rna_ids, hla_allele_type) %>% mutate(nb_tumour = n_distinct(tumour_id)) %>% mutate(proportion = nb_tumour/354) %>% group_by(symbol) %>% summarize(average = max(proportion), editing_level = mean(all_vaf)/2) %>% arrange(desc(average))

test$symbol <- factor(test$symbol, levels= as.character(test$symbol))

test %>% distinct() %>% pivot_longer(c(average, editing_level), names_to = "Type",values_to = "value") %>% ggplot(aes(x=symbol, y = value,fill= Type)) + 
  geom_col(position="dodge") + 
  scale_y_continuous(sec.axis = sec_axis(~ .*2, name = "Mean editing level"))+
  labs(y="Tumours with RNA editing (Proportion)", x= "Gene symbol") + my_theme + scale_x_discrete(guide = guide_axis(angle = 90)) + scale_fill_manual(name = "", labels = c("Tumours with RNA editing (Proportion)", "Mean editing level"), values = c("#F8766D","#00BFC4"))

test2 <- test %>% distinct() %>% pivot_longer(c(average, editing_level), names_to = "Type",values_to = "value")
test2$Type <- as.factor(test2$Type, levels = c("average","editing_level"))
test2 %>% ggplot(aes(x=value, y = reorder(symbol, value),fill= Type)) + 
  geom_col(position="dodge") + 
  scale_x_continuous(sec.axis = sec_axis(~ .*2, name = "Mean editing level")) +
  labs(y="Gene symbol", x= "Tumours with RNA editing (Proportion)") + my_theme + scale_y_discrete(guide = guide_axis(angle = 0)) + scale_fill_manual(name = "", labels = c("Tumours with RNA editing (Proportion)","Mean editing level"), values = c("#F8766D","#00BFC4"))

# Percentage of muttype -- Bar plot 
neoantigen <- merged_data %>% filter(effect == "nonsynonymous", EL_score < 0.1, EL_GL_score > EL_score, all_cov > 60)

neoantigen %>% group_by(muttype_nostrand) %>% summarize(no.rna_ids = n_distinct(rna_ids), proportion = n_distinct(rna_ids)*100/33) %>% arrange(desc(no.rna_ids)) %>% 
  ggplot(aes(x = muttype_nostrand,  y = proportion, fill = muttype_nostrand))+ geom_bar(stat = "identity") + 
  theme_bw() + labs(x = "Variant type", y = "Proportion of RNA editing events(%)") + 
  geom_text(aes(label = round(proportion, 1)), position = position_dodge(0.9),vjust = 1.5, color = "black") + my_theme2

# Compare the expression of ADAR and APOBEC3A between normal vs tumour regions ------------------------------------
neo_transcript3$sample_type <- as.factor(neo_transcript3$sample_type, levels = c("normal","tumour"))
neo_transcript3 %>% group_by(patient) %>% mutate(nb = n_distinct(sample_type)) %>% filter(nb > 1) %>% filter(gene_id== "ADAR") %>% group_by(patient, sample_type) %>% 
  summarize(expression = mean(expression)) %>% ggplot(aes(x = sample_type, y = expression, fill = sample_type)) + 
  geom_boxplot() + geom_point() + geom_line(aes(group= patient), linetype = "dashed") + 
  stat_compare_means() + theme_bw() + labs(x = "Sample type", y = "ADAR expression") + scale_fill_manual(values = normal_tumour_color)+ my_theme2

neo_transcript3 %>% group_by(patient) %>% mutate(nb = n_distinct(sample_type)) %>% filter(nb > 1) %>% filter(gene_id== "APOBEC3A") %>% group_by(patient, sample_type) %>% 
  summarize(expression = mean(expression)) %>% ggplot(aes(x = sample_type, y = expression, fill = sample_type)) + 
  geom_boxplot() + geom_point() + geom_line(aes(group= patient), linetype = "dashed") + 
  stat_compare_means() + theme_bw() + labs(x = "Sample type", y = "APOBEC3A expression") + scale_fill_manual(values = normal_tumour_color) + my_theme2

# Correlation between ADAR and neo-antigenic RNA editing variants ----------------------------------------------
neo_ids <- merged_data %>% filter(neoantigenic == TRUE)
neo_ids <- unique(neo_ids$rna_ids)

merged_data$neo <- FALSE 
merged_data$neo[which(merged_data$rna_ids %in% neo_ids)] <- TRUE

neo <- merged_data[which(merged_data$rna_ids %in% neo_ids),] %>% group_by(patient_region) %>% summarize(no_rna = n_distinct(rna_ids), no_neo = n_distinct(combine)) %>% column_to_rownames(var = "patient_region")
neo_merged <- merge(neo, ADAR_APO_exp, by = 0)
neo_merged %>% ggplot(aes(x = no_rna, y = ADAR)) + geom_point() + geom_smooth(method = "lm") + stat_cor(method = "spearman") + my_theme + labs(x = "No. of neo-antigenic variants (per region)", y = "ADAR expression")
neo_merged %>% ggplot(aes(x = no_rna, y = APOBEC3A)) + geom_point() + geom_smooth(method = "lm") + stat_cor(method = "spearman") + my_theme + labs(x = "No. of neo-antigenic variants (per region)", y = "APOBEC3A expression")

neo_proportion <- merged_data %>% group_by(patient_region) %>% mutate(total_rna = n_distinct(rna_ids)) %>% filter(neo == TRUE) %>% summarize(proportion = n_distinct(rna_ids)/total_rna) %>% distinct() %>% column_to_rownames(var = "patient_region")
neo_merged <- merge(neo_proportion, ADAR_APO_exp, by = 0)
neo_merged %>% ggplot(aes(x = proportion, y = ADAR)) + geom_point() + geom_smooth(method = "lm") + stat_cor(method = "spearman") + my_theme + labs(x = "Porportion of neo-antigenic variants", y = "ADAR expression")
neo_merged %>% ggplot(aes(x = proportion, y = APOBEC3A)) + geom_point() + geom_smooth(method = "lm") + stat_cor(method = "spearman") + my_theme + labs(x = "Porportion of neo-antigenic variants)", y = "APOBEC3A expression")

# Survival analysis ---------------------------------------------------------------------------------
## Import data 
survival <- readRDS("20220407_TRACERx421_all_tumour_df.rds")
colnames(survival)[1] <- "patient_region"
survival_data <- survival %>% select(patient_region, clinical_sex, patient_id, dfs_time, cens_dfs, os_time, cens_os) %>% distinct()
colnames(survival_data)[3] <- "patient"

## Count the no. of total RNA editing and no./proportion of neoantigenic RNA editing events per patient 
RE_burden <- merged_data %>% filter(!is.na(EL_score), effect =="nonsynonymous") %>% group_by(patient, patient_region) %>% mutate(total_rna = n_distinct(rna_ids), total_combine = n_distinct(combine)) %>% filter(neoantigenic == "TRUE") %>% mutate(neoantigen_burden = n_distinct(combine)) %>% select(patient_region,total_rna, total_combine, neoantigen_burden) %>% distinct()

## Merge with survival data 
merge_survival <- inner_join(RE_burden, survival_data)

# Cox proportional regression model 
survival_dfs <- Surv(merge_survival$dfs_time, merge_survival$cens_dfs)
survival_os <- Surv(merge_survival$os_time, merge_survival$cens_os)

coxph_total_rna <- coxph(survival_os ~ merge_survival$total_rna)
summary(coxph_total_rna)
coxph_total_combine <- coxph(survival_os ~ merge_survival$total_combine)
summary(coxph_total_combine)
coxph_neoantigen_burden <- coxph(survival_os ~ merge_survival$neoantigen_burden)
summary(coxph_neoantigen_burden)

# KM plot (RE neoantigen burden and OS)
merge_survival$variable_of_interest <- NA
merge_survival$variable_of_interest <- merge_survival$neoantigen_burden > median(merge_survival$neoantigen_burden,na.rm = T)

kmfit <- survfit(Surv(os_time, cens_os) ~ variable_of_interest, data= merge_survival)
ggsurvplot(kmfit, data = merge_survival, pval = TRUE, pval.method = TRUE, xlab= "Time (Days)",
           ylab= "Survival probability", title = "OS by RE neoantigen burden", ggtheme = theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()))

kmfit <- survfit(Surv(dfs_time, cens_dfs) ~ variable_of_interest, data= merge_survival)
ggsurvplot(kmfit, data = merge_survival, pval = TRUE, pval.method = TRUE, xlab= "Time (Days)",ylab= "Survival probability", 
            title = "DFS by RE neoantigen burden", ggtheme = theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()))


ggsurvplot(kmfit, data= merge_survival,
           pval = TRUE, pval.method = TRUE, conf.int = FALSE,
           risk.table = FALSE, # Add risk table
           risk.table.col = "strata", # Change risk table color by group 
           #surv.median.line = "hv", # Specify median survival
           legend.title = "",
           xlab= "Overall survival (Days)",
           ylab= "Survival probability",
           title = "OS by RE neoantigen burden",
           ggtheme = theme_bw() + my_theme)

# Correlate survival with ADAR & APOBEC3A expression
ADAR_APO_exp$patient_region <- rownames(ADAR_APO_exp)
merge_ADAR_survival <- inner_join(ADAR_APO_exp, survival_data)

# Cox proportional regression model 
survival_dfs <- Surv(merge_ADAR_survival$dfs_time, merge_ADAR_survival$cens_dfs)
survival_os <- Surv(merge_ADAR_survival$os_time, merge_ADAR_survival$cens_os)

coxph_ADAR <- coxph(survival_os ~ merge_ADAR_survival$ADAR)
summary(coxph_ADAR)
coxph_APOBEC3A <- coxph(survival_os ~ merge_ADAR_survival$APOBEC3A)
summary(coxph_APOBEC3A)

# KM plot (ADAR & APOBEC3A expression)
merge_ADAR_survival$variable_of_interest <- NA
merge_ADAR_survival$variable_of_interest <- merge_ADAR_survival$ADAR > median(merge_ADAR_survival$ADAR,na.rm = T)

kmfit <- survfit(Surv(os_time, cens_os) ~ variable_of_interest, data= merge_ADAR_survival)
ggsurvplot(kmfit, data = merge_ADAR_survival, pval = TRUE, pval.method = TRUE, xlab= "Time (Days)",
           ylab= "Survival probability", title = "OS by ADAR expression", ggtheme = theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()))
kmfit <- survfit(Surv(dfs_time, cens_dfs) ~ variable_of_interest, data= merge_survival)
ggsurvplot(kmfit, data = merge_survival, pval = TRUE, pval.method = TRUE, xlab= "Time (Days)",ylab= "Survival probability", 
           title = "DFS by ADAR expression", ggtheme = theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()))

merge_ADAR_survival$variable_of_interest <- NA
merge_ADAR_survival$variable_of_interest <- merge_ADAR_survival$APOBEC3A > median(merge_ADAR_survival$APOBEC3A,na.rm = T)
kmfit <- survfit(Surv(os_time, cens_os) ~ variable_of_interest, data= merge_ADAR_survival)
ggsurvplot(kmfit, data = merge_ADAR_survival, pval = TRUE, pval.method = TRUE, xlab= "Time (Days)",
           ylab= "Survival probability", risk.table = TRUE, title = "OS by APOBEC3A expression", ggtheme = theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()))
kmfit <- survfit(Surv(dfs_time, cens_dfs) ~ variable_of_interest, data= merge_survival)
ggsurvplot(kmfit, data = merge_survival, pval = TRUE, pval.method = TRUE, xlab= "Time (Days)",ylab= "Survival probability", 
           title = "DFS by APOBEC3A expression", ggtheme = theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()))

# Fig3 Compare between hot vs cold regions ------------------------------------------------
# Compare between hot vs cold regions -------------------------------------------------------------------------
## Mixed effect model-- no. total rna editing events & no. neoantigenic editing events 
immune_total_rna <- merged_data %>% filter(!is.na(immune_cat)) %>% group_by(immune_cat, tumour_id, patient, patient_region) %>% summarize(total_rna = n_distinct(rna_ids)) 
mixed.model_total.rna <- lmer(formula = total_rna ~ immune_cat + (1|tumour_id), data = immune_total_rna)
summary(mixed.model_total.rna)
  # Plot 
immune_total_rna %>% ggplot(aes(x = immune_cat, y = total_rna, fill =immune_cat)) + geom_boxplot()  + geom_point() + stat_compare_means() + theme_bw() + labs(x = "Immune category", y = "No. of total RNA editing variants") + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + theme(legend.position = "none") + my_theme2

immune_neo_rna <- merged_data %>% filter(!is.na(immune_cat), neoantigenic == TRUE) %>% group_by(immune_cat, tumour_id, patient, patient_region) %>% summarize(neo_rna = n_distinct(rna_ids)) 
mixed.model_neo.rna <- lmer(formula = neo_rna ~ immune_cat + (1|tumour_id), data = immune_neo_rna)
summary(mixed.model_neo.rna)
  # Plot 
immune_neo_rna %>% ggplot(aes(x = immune_cat, y = neo_rna, fill = immune_cat)) + geom_boxplot()  + geom_point() + stat_compare_means() + theme_bw() + labs(x = "Immune category", y = "Neo-antigenic RNA editing variants") + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + theme(legend.position = "none")+ my_theme2

# Mean editing level 
immune_editing_level <- merged_data %>% filter(!is.na(immune_cat), all_cov > 60, effect == "nonsynonymous") %>% group_by(immune_cat, tumour_id, patient, patient_region) %>% summarize(mean_editing_level = mean(all_vaf)) 
  #Plot 
immune_editing_level %>% ggplot(aes(x = immune_cat, y = mean_editing_level, fill =immune_cat)) + geom_boxplot()  + geom_point() + stat_compare_means() + theme_bw() + labs(x = "Immune category", y = "Mean editing level(VAF)") + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + theme(legend.position = "none") + my_theme2


# Compare no. nonsynonymous RNA ids, editing level & binding affinity between hot vs cold 
merged_data2 <- inner_join(merged_data, neoantigen2)
# No. nonsynonymous rna_ids 
merged_data2 %>% filter(!is.na(immune_cat), all_cov > 60, effect == "nonsynonymous") %>% group_by(immune_cat, patient, tumour_id, patient_region) %>% 
  summarize(total_rna = n_distinct(rna_ids), mean_editing_level = mean(all_vaf), average_binding_affinity = mean(BA_score)) %>% 
  ggplot(aes(x = immune_cat, y = total_rna, fill = immune_cat)) + 
  geom_boxplot()  + geom_point() + stat_compare_means() + theme_bw() + labs(x = "Immune category", y = "No. of nonsynonymous RNA editing variants") + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + theme(legend.position = "none")

# Average binding affinity 
merged_data2 %>% filter(!is.na(immune_cat), all_cov > 60, !is.na(EL_score), effect == "nonsynonymous") %>% group_by(immune_cat, patient, tumour_id, patient_region) %>% 
  summarize(mean_editing_level = mean(all_vaf), average_binding_affinity = mean(BA_score)) %>% 
  ggplot(aes(x = immune_cat, y = average_binding_affinity, fill = immune_cat)) + 
  geom_boxplot()  + geom_point() + stat_compare_means() + theme_bw() + labs(x = "Immune category", y = "Average BA score") + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + theme(legend.position = "none")

# Mean editing level 
merged_data %>% filter(!is.na(immune_cat), all_cov > 60, !is.na(EL_score)) %>% group_by(immune_cat, patient, tumour_id, patient_region) %>% 
  summarize(mean_editing_level = mean(all_vaf), average_binding_affinity = mean(BA_score)) %>% 
  ggplot(aes(x = immune_cat, y = mean_editing_level, fill = immune_cat)) + 
  geom_boxplot()  + geom_point() + stat_compare_means() + theme_bw() + labs(x = "Immune category", y = "Mean editing level") + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + theme(legend.position = "none")+ my_theme2

# No. of RNA ids by muttype_nostrand 
merged_data2 %>% filter(!is.na(immune_cat),muttype_nostrand != "C>A", muttype_nostrand != "C>G") %>% group_by(immune_cat,patient_region, muttype_nostrand) %>% 
  summarize(total_rna = n_distinct(rna_ids)) %>% ggplot(aes(x = muttype_nostrand, y = total_rna, fill = immune_cat)) + 
  geom_boxplot() + stat_compare_means(aes(label = ..p.signif..), label.y = 27) + theme_bw() + labs(x = "Variant type", y = "No. of RNA editing variants") + my_theme2

# Mixed effect model 
merged_data_mixed <- merged_data2 %>% filter(!is.na(immune_cat), all_cov > 60, !is.na(EL_score), effect == "nonsynonymous") %>% group_by(immune_cat, patient, tumour_id, patient_region) %>% 
  summarize(total_rna = n_distinct(rna_ids), mean_editing_level = mean(all_vaf), average_binding_affinity = mean(BA_score))

mixed_binding_affinity <- lmer(formula = average_binding_affinity ~ immune_cat + (1|tumour_id), data = merged_data_mixed)
summary(mixed_binding_affinity)

mixed_total_rna <- lmer(formula = total_rna ~ immune_cat + (1|tumour_id), data = merged_data_mixed)
summary(mixed_total_rna)

mixed_editing_level <- lmer(formula = mean_editing_level ~ immune_cat + (1|tumour_id), data = merged_data_mixed)
summary(mixed_editing_level)

# Compare the ADAR and APOBEC3A expression between hot vs cold regions 
neo_transcript3$immune_cat <- factor(neo_transcript3$immune_cat, levels = c("hot","cold"))
neo_transcript3 %>% filter(info != "neoantigens", gene_id == "ADAR", !is.na(immune_cat)) %>% # hot vs cold 
  ggplot(aes(x = immune_cat, y = expression, fill = immune_cat)) + geom_boxplot() + stat_compare_means() + theme_bw() + labs(y = "ADAR expression", x = "Immune category") + 
  geom_jitter(color = "grey", size = 0.7, alpha=0.5) + theme(legend.position = "none") + my_theme2

neo_transcript3 %>% filter(info != "neoantigens", gene_id == "APOBEC3A", !is.na(immune_cat)) %>% # hot vs cold 
  ggplot(aes(x = immune_cat, y = log2(expression), fill = immune_cat)) + geom_boxplot() + stat_compare_means() + theme_bw() + labs(y = "log2(APOBEC3A expression)") + 
  geom_jitter(color = "grey", size = 0.7, alpha=0.5) + theme(legend.position = "none") + my_theme2 + labs(x = "Immune category", y = "log2(APOBEC3A expression)")

# Compare NAGs genes expression between hot vs cold regions 
neo_transcript3 %>% filter(info == "neoantigens", !is.na(immune_cat)) %>% # hot vs cold 
  ggplot(aes(x = gene_id, y = log2(expression), fill = immune_cat)) + geom_boxplot() + stat_compare_means(aes(label = ..p.signif..), label.y = 10) + theme_bw() + 
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5),axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 14)) + labs(x = "Gene symbol")


# Correlation between immune cell abundance -------------------------------------------------------------------------
## NE neoantigen burden
neo_burden <- merged_data %>% group_by(patient_region) %>% mutate(total = n_distinct(combine)) %>% filter(neoantigenic == TRUE) %>% group_by(patient_region) %>% summarize(nb_neo = n_distinct(combine), neo_proportion = n_distinct(combine)/total) %>% distinct() %>% 
  column_to_rownames(var = "patient_region")

imm_abundance <- danaher_scores %>% select(patient_region, cell_type_title, score) %>% pivot_wider(names_from = cell_type_title, values_from = score) %>% distinct() %>% column_to_rownames(var = "patient_region")
imm_neo_merged <- merge(neo_burden, imm_abundance, by = 0)

imm_neo_merged %>% pivot_longer(c(4:19),names_to = "cell_type", values_to = "score") %>% ggplot(aes(x = nb_neo, y = score)) + geom_point() + geom_smooth(method = "lm") + 
  stat_cor(method = "spearman") + facet_wrap(~ cell_type) + theme_bw()

imm_neo_merged %>% pivot_longer(c(4:19),names_to = "cell_type", values_to = "score") %>% ggplot(aes(x = neo_proportion, y = score)) + geom_point() + geom_smooth(method = "lm") + 
  stat_cor(method = "spearman") + facet_wrap(~ cell_type) + theme_bw() + labs(x = "Proportion of neo-antigenic variants", y = "Danaher scores") + my_theme2

## Total RNA editing variants 
total_rna <- merged_data %>% group_by(patient_region) %>% summarize(total = n_distinct(rna_ids)) %>% column_to_rownames(var = "patient_region")
imm_total_rna_merged <- merge(total_rna, imm_abundance, by = 0)

imm_total_rna_merged %>% pivot_longer(c(3:18),names_to = "cell_type", values_to = "score") %>% ggplot(aes(x = total, y = score)) + geom_point() + geom_smooth(method = "lm") + 
  stat_cor(method = "spearman") + facet_wrap(~ cell_type) + theme_bw()

## ADAR expression & APOBEC3A expression 
ADAR_APO_exp2 <- ADAR_APO_exp %>% select(ADAR, APOBEC3A)
ADAR_imm_merge <- merge(ADAR_APO_exp2, imm_abundance, by = 0)

ADAR_imm_merge %>% pivot_longer(c(4:19),names_to = "cell_type", values_to = "score") %>% ggplot(aes(x = ADAR, y = score)) + geom_point() + geom_smooth(method = "lm") + 
  stat_cor(method = "spearman") + facet_wrap(~ cell_type) + labs(x = "ADAR expression", y = "Danaher scores") + theme_bw() + my_theme2

ADAR_imm_merge %>% pivot_longer(c(4:19),names_to = "cell_type", values_to = "score") %>% ggplot(aes(x = log2(APOBEC3A), y = score)) + geom_point() + geom_smooth(method = "lm") + 
  stat_cor(method = "spearman") + facet_wrap(~ cell_type) + labs(x = "log2(APOBEC3A expression)", y = "Danaher scores") + theme_bw()

# Figure 4: RNA editing and tumour immune evasion strategies -----------------------------------------------------
## Compare the expression of host genes between tumour and matched normal samples --------------------
neo_transcript3 %>% filter(info == "neoantigens") %>% # normal vs tumour 
  ggplot(aes(x = gene_id, y = log2(expression), fill = sample_type)) + geom_boxplot() + stat_compare_means(aes(label = ..p.signif..), label.y = 10) + theme_bw() + 
  scale_fill_manual(values = normal_tumour_color) + theme(axis.text.x = element_text(size = 9, angle = 120, vjust = 0.5)) + theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),axis.text.y = element_text(size = 10), 
                                                                                                                                  axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12)) + 
  labs(y ="log2(Expression)", x = "Gene symbol")

## Compare gene expression between patients with and without particular variants 
neo_ids <- unique(neoantigen$rna_ids)

neoantigen %>% group_by(symbol) %>% summarize(nb = n_distinct(rna_ids)) %>% filter(nb > 1 )
merged_data %>% filter(symbol == "CCNI"|symbol == "FLNA"|symbol == "FLNB") %>% group_by(symbol) %>% summarize(nb = n_distinct(rna_ids)) %>% filter(nb > 1 )

meta_data <- merged_data %>% ungroup() %>% select(patient_region) %>% distinct()

for (i in 1:length(neo_ids)){
  data <- merged_data %>% group_by(patient_region) %>% summarize(variant_present = ifelse(grepl(neo_ids[i], rna_ids),TRUE, FALSE)) %>% distinct()
  true <- data$patient_region[which(data$variant_present == TRUE)]
  data$variant_present[which(data$patient_region %in% true)] <- TRUE
  data <- unique(data)
  
  colnames(data)[2] <- unique(merged_data$symbol[which(merged_data$rna_ids %in% neo_ids[i])])
  meta_data <- left_join(meta_data, data)
}

# CCNI, FLNA and FLNB genes have two distinct RNA editing variants 
CCNI_true <- neoantigen %>% filter(symbol == "CCNI")
CCNI_true <- unique(CCNI_true$patient_region)

FLNA_true <- neoantigen %>% filter(symbol == "FLNA")
FLNA_true <- unique(FLNA_true$patient_region)

FLNB_true <- neoantigen %>% filter(symbol == "FLNB")
FLNB_true <- unique(FLNB_true$patient_region)

meta_data$CCNI <- FALSE
meta_data$CCNI[which(meta_data$patient_region %in% CCNI_true)] <- TRUE
meta_data$FLNA <- FALSE
meta_data$FLNA[which(meta_data$patient_region %in% FLNA_true)] <- TRUE
meta_data$FLNB <- FALSE
meta_data$FLNB[which(meta_data$patient_region %in% FLNB_true)] <- TRUE


# Merge with gene expression data 
meta_transcript_merge <- merge(neo_transcript3, meta_data)

gene_ids <- unique(neoantigen$symbol)
gene_ids <- gene_ids[-which(gene_ids == "HIST2H2BE")] # this gene is absent from transcript data 
gene_ids <- gene_ids[-which(gene_ids == "PLD6")] # this gene is absent from the meta data 

meta_transcript_merge2 <- data.frame()
for (i in 1:length(gene_ids)){
  data <- meta_transcript_merge %>% filter(gene_id == gene_ids[i]) %>% select(patient_region, gene_ids[i], gene_id, expression)
  colnames(data) <- c("patient_region","variant_present","gene_id","expression")
  data$group <- gene_ids[i]
  meta_transcript_merge2 <- rbind(meta_transcript_merge2, data)
}

a <- unique(meta_transcript_merge2$group)
a
meta_transcript_merge2$group <- as.factor(meta_transcript_merge2$group)
meta_transcript_merge2 %>% filter(group == a[1:12]) %>% group_by(group) %>% ggplot(aes(x = variant_present, y = log2(expression), fill = variant_present)) + geom_boxplot() + stat_compare_means(aes(label = ..p.signif..),label.x = 1.5) + theme_bw() + 
  theme(axis.text.x = element_text(size = 8, angle = 0, vjust = 0.5),axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 12), axis.title.x = element_blank()) + facet_wrap(~ group,scales = "free_y") + geom_jitter(color = "grey", size = 0.7, alpha=0.5)

meta_transcript_merge2 %>% filter(group == a[13:29]) %>% group_by(group) %>% ggplot(aes(x = variant_present, y = log2(expression), fill = variant_present)) + geom_boxplot() + stat_compare_means(aes(label = ..p.signif..), label.y = 6, label.x = 1.5) + theme_bw() + 
  theme(axis.text.x = element_text(size = 8, angle = 0, vjust = 0.5),axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 12), axis.title.x = element_blank()) + facet_wrap(~ group,scales = "free_y") + geom_jitter(color = "grey", size = 0.7, alpha=0.5)


meta_transcript_merge2 %>% group_by(group) %>% ggplot(aes(x = group, y = log2(expression), fill = variant_present)) + geom_boxplot() + stat_compare_means(aes(label = ..p.signif..),label.x = 1.5) + theme_bw() + 
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 12), axis.title.x = element_blank())  

## 
meta_transcript_merge2$variant_present <- factor(meta_transcript_merge2$variant_present, levels =c("FALSE","TRUE"))
meta_transcript_merge2 %>% filter(group == "LIMD1"|group == "FLNB"|group == "FLNA"|group == "COG3"|group == "EMP1"|
                                    group == "CTDNEP1"|group == "TGFBI"|group == "CCNI"|group == "EEF1G"|group == "PLEC"|
                                    group == "NCSTN"|group == "CDK12"|group == "NINL"|group == "TCF3"|group == "DCAF5"|
                                    group == "CYFIP2") %>% group_by(group) %>% ggplot(aes(x = variant_present, y = log2(expression), fill = variant_present)) + 
  geom_boxplot() + theme_bw() + facet_wrap(~ group, scales = "free_y") + 
  theme(axis.text.x = element_text(size = 8, angle = 0, vjust = 0.5),axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 12), axis.title.x = element_blank()) + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + geom_boxplot() + 
  scale_fill_manual(values = variant_present_color) + my_theme2

meta_transcript_merge2 %>% filter(group == "PARP10"|group == "TOR1B"|group == "TM7SF3"|group == "PABPC4"|group == "BLCAP"|
                                    group == "HES6"|group == "SLC35E3"|group == "CDK13"|group == "PAPSS1"|group == "KMT2D"|
                                    group == "FLCN"|group == "ODF2L"|group == "RPS6KC1") %>% group_by(group) %>% ggplot(aes(x = variant_present, y = log2(expression), fill = variant_present)) + 
  geom_boxplot() + theme_bw() + facet_wrap(~ group, scales = "free_y") + 
  theme(axis.text.x = element_text(size = 8, angle = 0, vjust = 0.5),axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 12), axis.title.x = element_blank()) + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + geom_boxplot() + 
  scale_fill_manual(values = variant_present_color) + my_theme2


meta_transcript_merge2 %>% group_by(group) %>% ggplot(aes(x = group, y = log2(expression), fill = variant_present)) + 
  geom_boxplot() + stat_compare_means(aes(label = ..p.signif..)) + theme_bw() + 
  theme(axis.text.x = element_text(size = 8, angle = 0, vjust = 0.5),axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 12), axis.title.x = element_blank()) 


meta_transcript_merge2 %>% filter(group == "PLEC") %>% distinct() %>% ggplot(aes(x= variant_present, y = log2(expression), fill = variant_present)) + geom_boxplot() + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + 
  stat_compare_means(aes(label = ..p.signif..), label.y = 8, label.x = 1.5) + theme_bw() + my_theme2 + labs(title = "PLEC") + 
  theme(axis.text.x = element_text(size = 8, angle = 0, vjust = 0.5),axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 12), axis.title.x = element_blank())


# Compare the expression of genes with NAGs between hot vs cold regions & between normal vs tumour regions 
neo_transcript3 %>% filter(info == "neoantigens", !is.na(immune_cat)) %>% # hot vs cold 
  ggplot(aes(x = gene_id, y = log2(expression), fill = immune_cat)) + geom_boxplot() + stat_compare_means(aes(label = ..p.signif..), label.y = 10) + theme_bw() + 
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5),axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 14))

## Compare the HLA LOH pattern between patients with and without HLA LOH ---------------------------------------
load("merged_lohhla.RData") 

# LOHHLA by patient level 
# LOHHLA true patients (n= 117), LOHHLA false patients (n= 210)
lohhla_patient <- merged_lohhla %>% filter(loh == TRUE)
lohhla_patient <- unique(lohhla_patient$patient)
merged_lohhla$LOH_by_patient <- FALSE
merged_lohhla$LOH_by_patient[which(merged_lohhla$patient %in% lohhla_patient)] <- TRUE

# Comparison of RNA editing pattern between patients with and without HLA LOH 
merged_lohhla$LOH_by_patient <- factor(merged_lohhla$LOH_by_patient, levels = c(TRUE, FALSE))

total.rna <- merged_lohhla %>% filter(effect == "nonsynonymous") %>%  group_by(patient, LOH_by_patient) %>% summarize(frequency = n_distinct(combine),mean_editing_level = mean(all_vaf)) %>% ggplot(aes(x = LOH_by_patient, y = frequency, fill = LOH_by_patient))  + stat_compare_means() + geom_boxplot() + geom_point() + labs(x = "LOHHLA", y = "Non-synonymous RNA editing variants") + theme_bw() + scale_fill_manual(values = lohhla_color) + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + theme(legend.position = "none") + my_theme2

neo.rna <- merged_lohhla %>% filter(all_cov > 60, effect == "nonsynonymous",EL_score < 0.1, EL_score < EL_GL_score) %>% group_by(patient, LOH_by_patient) %>% summarize(neoantigen = n_distinct(combine),mean_editing_level = mean(all_vaf)) %>% ggplot(aes(x = LOH_by_patient, y = neoantigen, fill = LOH_by_patient))  + stat_compare_means() + geom_boxplot()  + labs(x = "LOHHLA", y = "RE neoantigen burden") + theme_bw() + scale_fill_manual(values = lohhla_color) + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + theme(legend.position = "none")+ my_theme2

level <- merged_lohhla %>% filter(all_cov > 60, effect == "nonsynonymous") %>% group_by(patient, LOH_by_patient) %>% summarize(frequency = n_distinct(rna_ids),mean_editing_level = mean(all_vaf)) %>% ggplot(aes(x = LOH_by_patient, y = mean_editing_level, fill = LOH_by_patient))  + stat_compare_means() + geom_boxplot()  + labs(x = "LOHHLA", y = "Mean editing level") + theme_bw()+ scale_fill_manual(values = lohhla_color) + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + theme(legend.position = "none")+ my_theme2

EL_score <- merged_lohhla %>% filter(all_cov > 60, effect == "nonsynonymous") %>% group_by(patient, LOH_by_patient) %>% summarize(Average_EL = mean(EL_score)) %>%  ggplot(aes(x = LOH_by_patient, y = Average_EL, fill= LOH_by_patient)) + geom_boxplot() + geom_point() + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + stat_compare_means() + labs(x = "LOHHLA", y = "Average EL score") + theme_bw() + theme(legend.position = "none") + scale_fill_manual(values = lohhla_color) + my_theme2

BA_score <- merged_lohhla %>% filter(all_cov > 60, effect == "nonsynonymous") %>% group_by(patient, LOH_by_patient) %>% summarize(Average_BA = mean(BA_score)) %>%  ggplot(aes(x = LOH_by_patient, y = Average_BA, fill= LOH_by_patient)) + geom_boxplot() + geom_point() + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + stat_compare_means() + labs(x = "LOHHLA", y = "Average BA score") + theme_bw() + theme(legend.position = "none") + scale_fill_manual(values = lohhla_color) + my_theme2

neo_proportion <- merged_lohhla %>% group_by(patient) %>% mutate(total = n_distinct(rna_ids)) %>% filter(effect == "nonsynonymous", EL_score < 0.1, all_cov > 60 ) %>% group_by(patient, LOH_by_patient) %>% summarize(proportion = n_distinct(rna_ids)/total) %>% distinct() %>% ggplot(aes(x = LOH_by_patient, y = proportion, fill = LOH_by_patient)) + stat_compare_means() + geom_boxplot() + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + labs(x = "LOHHLA", y = "Proportion of neoantigenic RNA editing events", title = "Comparison of RNA editing proportion") + theme_bw() + theme(legend.position = "none") + scale_fill_manual(values = lohhla_color)

ggarrange(proportion, level, labels = c("A", "B"),ncol = 2, nrow = 1)

# Identify which HLA allele type was lost 
allele1 <- lohhla_data %>% filter(loh == TRUE, allele_type == "allele.1")
allele1$HLA_loss <- ifelse(allele1$cn1_binned < allele1$cn2_binned, "TRUE","FALSE")

allele2 <- lohhla_data %>% filter(loh == TRUE, allele_type == "allele.2")
allele2$HLA_loss <- ifelse(allele2$cn2_binned < allele2$cn1_binned, "TRUE", "FALSE")

lohhla_true <- rbind(allele1, allele2) %>% arrange(patient_region, patient, allele1, allele2)
lohhla_true <- lohhla_true %>% select(patient_region, patient, tumour_id, hla_type, loh, allele_type, HLA, HLA_loss) %>% distinct()

# Get Neoantigenic RNA editing events 
neoantigen <- merged_data2 %>% filter(effect == "nonsynonymous", EL_score < 0.1, EL_GL_score > EL_score, all_cov > 60)
nonsyn <- merged_data2 %>% filter(effect == "nonsynonymous", EL_score < 10, all_cov > 60) %>% select(rna_ids, HLA, hla_type, hla_allele_type,combine, patient_region, patient,tumour_id, EL_score, BA_score, all_vaf) %>% distinct()

test <- inner_join(nonsyn, lohhla_true)

# Compare the binding affinity (EL_scores/BA_score) between neoantigens with HLA_loss and not loss 
test$HLA_loss <- factor(test$HLA_loss, levels = c("TRUE","FALSE"))
test %>% group_by(patient, HLA_loss) %>% summarize(average_EL_score = mean(EL_score), mean_editing = mean(all_vaf), average_binding = mean(BA_score)) %>% 
  ggplot(aes(x = HLA_loss, y = average_EL_score, fill = HLA_loss)) + geom_boxplot() + geom_point() + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + stat_compare_means() + theme_bw() + theme(legend.position = "none") + scale_fill_manual(values = lohhla_color)

test %>% group_by(patient, HLA_loss) %>% summarize(average_EL_score = mean(EL_score), mean_editing = mean(all_vaf), average_binding = mean(BA_score)) %>% 
  ggplot(aes(x = HLA_loss, y = mean_editing, fill = HLA_loss)) + geom_boxplot() + geom_point() + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + stat_compare_means() + theme_bw() + theme(legend.position = "none") + scale_fill_manual(values = lohhla_color)

test %>% group_by(patient, HLA_loss) %>% summarize(average_EL_score = mean(EL_score), mean_editing = mean(all_vaf), average_binding = mean(BA_score)) %>% 
  ggplot(aes(x = HLA_loss, y = average_binding, fill = HLA_loss)) + geom_boxplot() + geom_point() + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + stat_compare_means() + theme_bw() + theme(legend.position = "none") + scale_fill_manual(values = lohhla_color)

# Among patients with both TRUE and FALSE
test %>% group_by(patient) %>% mutate(nb = n_distinct(HLA_loss)) %>% filter(nb > 1) %>% 
  group_by(patient, HLA_loss) %>% summarize(average_EL = mean(EL_score)) %>% distinct() %>% 
  ggplot(aes(x = HLA_loss, y = average_EL, fill= HLA_loss)) + geom_boxplot() + geom_point() + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + stat_compare_means() + labs(x = "HLA loss", y = "EL score") + theme_bw() + theme(legend.position = "none") + scale_fill_manual(values = lohhla_color)+ my_theme2

test %>% group_by(patient) %>% mutate(nb = n_distinct(HLA_loss)) %>% filter(nb > 1) %>% ggplot(aes(x = HLA_loss, y = EL_score, fill= HLA_loss)) + geom_boxplot() + geom_point() + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + stat_compare_means() + labs(x = "HLA loss", y = "Binding affinity") + theme_bw() + theme(legend.position = "none") + scale_fill_manual(values = lohhla_color)

test %>% group_by(patient) %>% mutate(nb = n_distinct(HLA_loss)) %>% filter(nb > 1) %>% ggplot(aes(x = HLA_loss, y = BA_score, fill= HLA_loss)) + geom_boxplot() + geom_point() + geom_jitter(color = "grey", size = 0.7, alpha=0.5) + stat_compare_means() + labs(x = "HLA loss", y = "BA score") + theme_bw() + theme(legend.position = "none") + scale_fill_manual(values = lohhla_color)

test %>% group_by(patient) %>% mutate(nb = n_distinct(HLA_loss)) %>% filter(nb > 1) %>% ggplot(aes(x = HLA_loss, y = all_vaf)) + geom_boxplot() + geom_point() + stat_compare_means()

# Pie chart + Bar plot 
data <- data.frame(
  HLA_loss = c("TRUE", "FALSE"),
  No.patients = c(61, 135))

data <- data %>% mutate(perc = No.patients/sum(No.patients)) %>% arrange(perc) %>% mutate(labels = scales::percent(perc)) 

data2 <- data %>% mutate(csum = rev(cumsum(rev(perc*100))), 
                         pos = perc*100/2 + lead(csum, 1),
                         pos = if_else(is.na(pos), perc*100/2, pos))

ggplot(data, aes(x = "" , y = perc*100, fill = fct_inorder(HLA_loss))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  geom_label_repel(data = data2,
                   aes(y = pos, label = labels ),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "HLA LOH")) +
  theme_void() + scale_fill_manual(values = lohhla_color) 

# Plot the percentage of HLA loss by patient 
test %>% group_by(patient) %>% mutate(total_combine = n_distinct(combine)) %>% filter(HLA_loss == "TRUE") %>% group_by(patient) %>% mutate(proportion = n_distinct(combine)/total_combine) %>% select(patient, total_combine, proportion) %>%  distinct() %>% 
  ggplot(aes(x = proportion, y = reorder(patient, proportion))) + geom_bar(stat ="identity") + theme_bw() + my_theme2 + labs(x = "Proportion of disrupted RE neoantigen", y = "Patient")
  
  
  #geom_text(aes(label = paste0("n=",total_combine)), position = position_dodge(-1),
  #          vjust = 1.5, color = "black") 

# Survival analysis for patients with and without LOHHLA


merged_data %>% group_by(patient) %>% mutate(nb = n_distinct(sample_type)) %>% filter(nb > 1, all_cov > 60) %>% 
  group_by(patient, patient_region,sample_type) %>% 
  summarize(total_rna = n_distinct(rna_ids), mean_editing = mean(all_vaf)) %>% group_by(patient, sample_type) %>% 
  summarize(av_total_rna = mean(total_rna), av_mean_editing = mean(mean_editing)) %>% distinct() %>% 
  ggplot(aes(x = sample_type, y = av_total_rna, fill = sample_type )) + geom_boxplot() + geom_point() + stat_compare_means() + geom_line(aes(group = patient), linetype = "dashed")

merged_data %>% group_by(patient) %>% mutate(nb = n_distinct(sample_type)) %>% filter(nb > 1, all_cov > 60) %>% 
  group_by(patient, patient_region,sample_type) %>% 
  summarize(total_rna = n_distinct(rna_ids), mean_editing = mean(all_vaf)) %>% group_by(patient, sample_type) %>% 
  summarize(av_total_rna = mean(total_rna), av_mean_editing = mean(mean_editing)) %>% distinct() %>% 
  ggplot(aes(x = sample_type, y = av_mean_editing, fill = sample_type )) + geom_boxplot() + geom_point() + stat_compare_means() + geom_line(aes(group = patient), linetype = "dashed")


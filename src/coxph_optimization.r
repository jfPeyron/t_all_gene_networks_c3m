library(readxl)
library(survival)
library(survminer)
library(glmnet)
library(MASS)

setwd("/src")

# Global var -------------------------------------------------------------------
PRINT_GRP_TABLE = FALSE
PLOT_KAPLAN_MEIER = FALSE

# Set up cohort ----------------------------------------------------------------
expression <- as.data.frame(readRDS("../data/cohort/data.table.rds"))
expression <- as.data.frame(expression)
expression$`TARGET USI` <- rownames(expression)

patient_data <- as.data.frame(read_excel("../data/cohort/clinical_data.xlsx"))
patient_data <- patient_data[patient_data$Protocol == "AALL0434",]

patient_data.filtered <- patient_data[names(patient_data) %in% c("TARGET USI","Event Free Survival Time in Days","First Event")]
rownames(patient_data.filtered) <- rownames(patient_data)
patient_data.filtered.noNA <- na.omit(patient_data.filtered)
cohort_data <- merge(expression,patient_data.filtered.noNA, by="TARGET USI")


# Define vars extracted from sPLS-DA from each library -------------------------
lib1.vars <- readRDS("../data/cox_result/lib1_outsplsda.rds")
lib2.vars <- readRDS("../data/cox_result/lib2_outsplsda.rds")
lib3.vars <- readRDS("../data/cox_result/lib3_outsplsda.rds")
libs <- list(lib1.vars, lib2.vars, lib3.vars)


lib_id = 1
for (lib_vars in libs) {
  
  cohort_data.genes_expression.scoreapplied <- cohort_data[names(cohort_data) %in% lib_vars]
  cohort_data.genes_expression.scoreapplied <- cbind(cohort_data.genes_expression.scoreapplied,"OS" = cohort_data$`Event Free Survival Time in Days`)
  cohort_data.genes_expression.scoreapplied <- cbind(cohort_data.genes_expression.scoreapplied,"Vital" = cohort_data$`First Event`)
  
  cohort_data.genes_expression.scoreapplied$Vital <- gsub("None","0",cohort_data.genes_expression.scoreapplied$Vital)
  cohort_data.genes_expression.scoreapplied$Vital <- gsub("Relapse","1",cohort_data.genes_expression.scoreapplied$Vital)
  cohort_data.genes_expression.scoreapplied$Vital <- as.numeric(cohort_data.genes_expression.scoreapplied$Vital)
  
  fit <- coxph(Surv(`OS`,`Vital`) ~ .,data = cohort_data.genes_expression.scoreapplied)
  
  fit.select <- stepAIC(fit, direction = "both", trace = 1)
  
  coef_table <- summary(fit.select)$coefficients[, c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")]
  model_selection <- data.frame(coef_table)
  
  ## Define risk groups --------------------------------------------------------
  cohort_data.genes_expression.genesselected <- cohort_data.genes_expression.scoreapplied[names(cohort_data.genes_expression.scoreapplied) %in% names(coef(fit.select))]
  temps <- cohort_data.genes_expression.scoreapplied$OS
  censure <- cohort_data.genes_expression.scoreapplied$Vital
  score <- sweep(cohort_data.genes_expression.genesselected, 2, coef(fit.select), FUN = "*")
  score <- cbind(score, "score" = apply(score, 1, sum))
  groupe <- ifelse(score$score >= median(score$score), "High", "Low")
  donnees <- data.frame(temps, censure, groupe)
  
  if (PRINT_GRP_TABLE) {
    table_censure_groupe <- table(donnees$censure, donnees$groupe)
    
    print(paste("###### LIB ", lib_id, " ######", sep = ""))
    print(table_censure_groupe)
  }
  
  if (PLOT_KAPLAN_MEIER) {
    courbe_survie <- survfit(Surv(temps, censure) ~ groupe, data = donnees)
    p = ggsurvplot(courbe_survie, data = donnees, risk.table = T, pval = TRUE, conf.int = FALSE, break.y.by = 0.10)
    plot(p$plot)
  }
  
  lib_id <- lib_id + 1
}
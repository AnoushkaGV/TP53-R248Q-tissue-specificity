# Load required libraries
library(dplyr)
library(survival)
library(gridExtra)
library(grid)

# Read mutation and clinical data for multiple cancer types
brca_mut <- read.delim("./data/mutations_brca.txt")
brca_clin <- read.delim("./data/brca_tcga_pan_can_atlas_2018_clinical_data.tsv")

sarc_mut <- read.delim("./data/mutations_sarc.txt")
sarc_clin <- read.delim("./data/sarc_tcga_pan_can_atlas_2018_clinical_data.tsv")

gbm_mut <- read.delim("./data/mutations_gbm.txt")
gbm_clin <- read.delim("./data/gbm_tcga_pan_can_atlas_2018_clinical_data.tsv")

coad_mut <- read.delim("./data/mutations_coad.txt")
coad_clin <- read.delim("./data/coadread_tcga_pan_can_atlas_2018_clinical_data.tsv")

# Function to process mutation and clinical data
process_data <- function(mut, clin, cancer_name) {
  # Identify TP53 R248 mutations vs other TP53 mutations
  tp53 <- mut %>%
    mutate(patient = SAMPLE_ID,
           group = ifelse(grepl("R248", TP53), "R248", "Other_TP53")) %>%
    select(patient, group) %>%
    distinct()
  
  # Merge with clinical data and define survival variables
  clin %>%
    mutate(patient = Sample.ID) %>%
    left_join(tp53, by = "patient") %>%
    mutate(
      group = ifelse(is.na(group), "WT", group),   # Wildtype if no TP53 mutation
      time = as.numeric(Overall.Survival..Months.),  # Survival time in months
      status = ifelse(Overall.Survival.Status == "1:DECEASED", 1, 0), # Event indicator
      age = as.numeric(Diagnosis.Age),
      sex = Sex,
      cancer = cancer_name
    ) %>%
    filter(!is.na(time), !is.na(status), time >= 1) %>%
    select(patient, group, time, status, age, sex, cancer)
}

# Process each cancer type
brca <- process_data(brca_mut, brca_clin, "Breast")
sarc <- process_data(sarc_mut, sarc_clin, "Sarcoma")
gbm <- process_data(gbm_mut, gbm_clin, "Glioblastoma")
coad <- process_data(coad_mut, coad_clin, "Colorectal")

# Combine all cancer types into a single dataset
data <- bind_rows(brca, sarc, gbm, coad)
data$group <- factor(data$group, levels = c("WT", "Other_TP53", "R248"))

# ------------------------------
# Descriptive statistics
# ------------------------------
summary_overall <- data %>%
  group_by(cancer, group) %>%
  summarise(
    n = n(),
    median_age = round(median(age, na.rm = TRUE), 1),
    n_female = sum(sex == "Female", na.rm = TRUE),
    n_male = sum(sex == "Male", na.rm = TRUE),
    deaths = sum(status),
    median_OS = round(median(time), 1),
    .groups = "drop"
  )
print(summary_overall)

# Display descriptive table as a figure
grid.newpage()
grid.table(summary_overall, rows = NULL,
           theme = ttheme_default(base_size = 14, padding = unit(c(8,4), "mm")))

# Compare ages between R248 and WT patients
for (cancer in unique(data$cancer)) {
  d <- filter(data, cancer == !!cancer, !is.na(age))
  wt <- d$age[d$group == "WT"]
  r248 <- d$age[d$group == "R248"]
  if (length(r248) >= 5) {
    t_test <- t.test(r248, wt)
    cat(cancer, "- R248 vs WT: p =", format.pval(t_test$p.value, digits = 3), "\n")
  }
}

# ------------------------------
# Kaplan-Meier survival analysis
# ------------------------------
layout(matrix(c(1,2,3,4,5,5), nrow = 3, byrow = TRUE), heights = c(1,1,0.2))
par(mar = c(5,5,4,2), cex.lab = 1.4, cex.axis = 1.2, cex.main = 1.5)

for (cancer in unique(data$cancer)) {
  d <- filter(data, cancer == !!cancer)
  fit <- survfit(Surv(time, status) ~ group, data = d)
  
  logrank <- survdiff(Surv(time, status) ~ group, data = d)
  pval <- pchisq(logrank$chisq, length(logrank$n)-1, lower.tail = FALSE)
  
  cat(cancer, ": p =", format.pval(pval, digits = 3), "\n")
  
  plot(fit, col = c("#00BA38","#619CFF","#F8766D"), lwd = 3,
       xlab = "Time (months)", ylab = "Survival Probability",
       main = paste(cancer, "\np =", format.pval(pval, digits = 3)))
}

# Add legend at bottom
par(mar = c(0,0,0,0))
plot.new()
legend("center", legend = c("WT","Other TP53","R248"),
       col = c("#00BA38","#619CFF","#F8766D"), lwd = 4,
       bty = "n", cex = 1.5, ncol = 3, horiz = TRUE)

par(mfrow = c(1,1))
layout(1)

# ------------------------------
# Cox proportional hazards models
# ------------------------------
cox_results <- data.frame()
for (cancer in unique(data$cancer)) {
  d <- filter(data, cancer == !!cancer, !is.na(age))
  if (nrow(d) < 30) next
  
  model <- coxph(Surv(time, status) ~ group + age, data = d)
  s <- summary(model)
  
  if ("groupR248" %in% rownames(s$coefficients)) {
    cox_results <- rbind(cox_results, data.frame(
      cancer = cancer,
      comparison = "R248 vs WT",
      HR = s$conf.int["groupR248",1],
      lower = s$conf.int["groupR248",3],
      upper = s$conf.int["groupR248",4],
      p_value = s$coefficients["groupR248",5]
    ))
  }
  
  if ("groupOther_TP53" %in% rownames(s$coefficients)) {
    cox_results <- rbind(cox_results, data.frame(
      cancer = cancer,
      comparison = "Other TP53 vs WT",
      HR = s$conf.int["groupOther_TP53",1],
      lower = s$conf.int["groupOther_TP53",3],
      upper = s$conf.int["groupOther_TP53",4],
      p_value = s$coefficients["groupOther_TP53",5]
    ))
  }
}

cox_results$p_adj <- p.adjust(cox_results$p_value, method = "BH")
cox_results$significant <- ifelse(cox_results$p_adj < 0.05, "*", "")
print(cox_results)

# Display Cox regression results as table
grid.newpage()
cox_display <- cox_results %>%
  mutate(HR = round(HR,2),
         lower = round(lower,2),
         upper = round(upper,2),
         p_value = format.pval(p_value, digits = 2),
         p_adj = format.pval(p_adj, digits = 2)) %>%
  select(cancer, comparison, HR, lower, upper, p_value, p_adj, significant)
grid.table(cox_display, rows = NULL,
           theme = ttheme_default(base_size = 14, padding = unit(c(8,4), "mm")))

# ------------------------------
# Forest plot for R248 vs WT hazard ratios
# ------------------------------
forest_data <- cox_results %>% filter(comparison == "R248 vs WT")
par(mar = c(5,8,4,2), cex.lab = 1.4, cex.axis = 1.3, cex.main = 1.6)
plot(forest_data$HR, 1:nrow(forest_data), xlim = c(min(forest_data$lower)*0.8, max(forest_data$upper)*1.2),
     yaxt = "n", ylab = "", xlab = "Hazard Ratio", log = "x",
     pch = 19, col = "#F8766D", cex = 2.5,
     main = "R248 vs WT Hazard Ratios (adjusted for age)")
abline(v = 1, lty = 2, lwd = 2)
segments(forest_data$lower, 1:nrow(forest_data), forest_data$upper, 1:nrow(forest_data), lwd = 3)
axis(2, at = 1:nrow(forest_data), labels = forest_data$cancer, las = 1, cex.axis = 1.3)
text(forest_data$HR, 1:nrow(forest_data), 
     labels = paste0("HR=", round(forest_data$HR,2), forest_data$significant),
     pos = 3, cex = 1.2)


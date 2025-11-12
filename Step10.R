#STEP 10 
#11-03-2025

#reading in MAF files
brca <- read.delim("mutations_brca.txt", header = T, sep = "\t")
coad <- read.delim("mutations_coad.txt", header = T, sep = "\t")
gbm <- read.delim("mutations_gbm.txt", header = T, sep = "\t")
sarc <- read.delim("mutations_sarc.txt", header = T, sep = "\t")

#number of samples in each dataset
n_brca <- nrow(brca)
n_coad <- nrow(coad)
n_gbm <- nrow(gbm)
n_sarc <- nrow(sarc)

#number of unique variants in each dataset
nvar_brca <- length(unique(brca$TP53))
nvar_coad <- length(unique(coad$TP53))
nvar_gbm <- length(unique(gbm$TP53))
nvar_sarc <- length(unique(sarc$TP53))

#plotting mutation counts
#plotting frequency counts
library(ggplot2)

mut_plot <- function (df, sample_name){
  table1 <- table(df$TP53)
  table1 <- sort(table1, decreasing = T)
  
  mut_df <- as.data.frame(table1)
  colnames(mut_df) <- c("Mutation", "Count")
  
  mut_df1 <- subset(mut_df, Mutation != "WT")
  
  ggplot2::ggplot(head(mut_df1,20), aes(x = reorder(Mutation, - Count), y = Count)) +
    geom_bar(stat = "identity", fill = "pink") +
    theme_minimal() +
    labs(title = paste(sample_name),
         x = "Mutation", y = "Number of Samples") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

mut_plot(brca, "Top 20 TP53 Mutations in Invasive Breast Cancer Cohort")
mut_plot(coad, "Top 20 TP53 Mutations in Colorectal Adenocarcinoma Cohort")
mut_plot(gbm, "Top 20 TP53 Mutations in Glioblastoma Cohort")
mut_plot(sarc, "Top 20 TP53 Mutations in Sarcoma Cohort")

#isolate R248Q variant
r248_mut <- function(df, sample_name){
  #number of r248q mutations in the cohort
  r248q_samp <- which(df$TP53 == "R248Q")
  r248q <- df[r248q_samp,]
  nr248q <- nrow(r248q)
  cat(sprintf("The number of R248Q in the cohort: %d", nr248q))
  
  #other r248 mutations
  r248_samp <- grep("R248", df$TP53, value = F)
  r248 <- df[r248_samp, ]
  table1 <- as.data.frame(table(r248$TP53))
  colnames(table1) <- c("Mutation", "Count")
  table1 <- table1[order(-table1$Count), ]
  
  #frequency plot the R248 mutations
  ggplot2::ggplot(table1, aes(x = reorder(Mutation, -Count), y = Count)) +
    geom_bar(stat = "identity", fill = "hotpink") +
    labs(title = sprintf("Distribution of R248 mutations in %s", sample_name),
         x = "Mutation", y = "Number of Samples") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

r248_mut(brca, "Invasive Breast Cancer Cohort")
r248_mut(coad, "Colorectal Adenocarcinoma Cohort")
r248_mut(gbm, "Gliobastoma Cohort")
r248_mut(sarc, "Sarcoma Cohort")

#distribution of R248Q across different cancers
x <- c(5, 15, 6, 1)
y <- c("BRCA", "COAD", "GBM", "SARC")
mut1 <- data.frame(y,x)
head(mut1)
colnames(mut1) <- c("Cancer_type", "Count")

ggplot2::ggplot(mut1, aes(x= reorder(Cancer_type, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "salmon") +
  labs(title = "Frequency of R248Q mutation across different cancer types",
       x = "Cancer type", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#distribution of number of variants across the samples
a <- c(nvar_brca, nvar_coad, nvar_gbm, nvar_sarc)
b <- c("BRCA", "COAD", "GBM", "SARC")

mut2 <- data.frame(b,a)
colnames(mut2) <- c("Cancer_type", "Count")
head(mut2)

ggplot2::ggplot(mut2, aes(x= reorder(Cancer_type, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  labs(title = "Distribution of number of TP53 variants across the cohorts",
       x = "Cancer type", y = "Number of variants") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


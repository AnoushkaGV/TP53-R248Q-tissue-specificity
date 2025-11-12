#NGS Analysis

#reading mutation file
mut <- read.delim("mutations.txt", header = T, sep = "\t")
head(mut)

#number of samples in the dataset
nsamp <- nrow(mut)
nsamp
#996 samples

#unique variants
nvar <- length(unique(mut$TP53))
nvar
#199 unique variants

#mutation frequency counts
par(mfrow = c(2,1))
table1 <- table(mut$TP53)
table1 <- sort(table1, decreasing = T)
View(table1)
plot(table1, 
     main = "Mutation frequency", 
     xlab = "Variants",
     ylab = "Frequency",
     col = "pink")

head(table1, 10)
plot(table1[1:11],
     main = "Top 10 variants",
     xlab = "Variants",
     ylab = "Frequency",
     col = "hotpink")

#plotting frequency counts
library(ggplot2)

mut_df <- as.data.frame(table1)
colnames(mut_df) <- c("Mutation", "Count")

mut_df1 <- subset(mut_df, Mutation != "WT")
mut_df$Count <- as.numeric(mut_df$Count)
par(mfrow = c(2,1))
#plot mutations
ggplot(mut_df1, aes(x = reorder(Mutation, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "pink") +
  theme_minimal() +
  labs(title = "Mutation frequency in TCGA Breast Cancer Cohort",
       x = "Mutation", y = "Number of Samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot top 10 mutations
ggplot(head(mut_df1, 10), aes(x = reorder(Mutation, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "pink") +
  theme_minimal() +
  labs(title = "Top 10 TP53 Mutations in TCGA Breast Cancer Cohort",
       x = "Mutation", y = "Number of Samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#isolate wild type
wtsamp <- which(mut$TP53 == "WT")
wt <- mut[wtsamp,]
head(wt)
nrow(wt)

#isolate R248Q variant
r248samp <- which(mut$TP53 == "R248Q")
r248q <- mut[r248samp,]
nr248q <- nrow(r248q);nr248q
head(r248q)

#other R248 mutations
r248_mutsamp <- grep("R248", mut$TP53, value = FALSE)
r248 <- mut[r248_mutsamp, ]
nrow(r248)
r248c <- as.data.frame(table(r248$TP53))
colnames(r248c) <- c("Mutation", "Count")
r248c <- r248c[order(-r248c$Count), ]
head(r248c)


#frequency of other r248 mutations
ggplot(r248c, aes(x = reorder(Mutation, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "hotpink") +
  labs(title = "Distribution of R248 mutations in TCGA Breast Cancer Cohort",
       x = "Mutation", y = "Number of Samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#mutation proportion
wt_count <- sum(mut$TP53 == "WT")
mut_count <- nrow(mut) - wt_count
prop <- mut_count/wt_count
prop

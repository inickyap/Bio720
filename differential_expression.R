library(DESeq2)
library("RColorBrewer")
library("gplots")

setwd("~/Desktop/counts/data/")

#Importing data
in_dir = dir(, pattern= "_htseq_counts.txt" )
counts_in <- lapply(in_dir, function(x) read.table(x, header=F, nrows=23337))
tot_count_matrix <- matrix(unlist(lapply(counts_in, function(x) x$V2)) , ncol=48, nrow=23337)
parse_names <- strsplit(in_dir, split="_")
parse_names <- matrix(unlist(parse_names), nrow=48, ncol=8, byrow=T)
col_names_counts <- paste(parse_names[,1], "_", parse_names[,2], "_", parse_names[,3], parse_names[,4], sep="")
colnames(tot_count_matrix) = col_names_counts 
rownames(tot_count_matrix) = counts_in[[1]]$V1

#Setting up experimental design 
experimental_design = data.frame(
  sample_names = col_names_counts,  # sample name
  age = factor(parse_names[,1]), # old or young 
  treatment = factor(parse_names[,2]), # treatment plan 
  lane = factor(parse_names[,4])      # Which lane on the Illumina flowcell.
)


#Testing for lane effect 
test_lane <- DESeqDataSetFromMatrix(tot_count_matrix, experimental_design, design = formula(~ lane))
test_lane <- DESeq(test_lane) 
test_lane_results <- results(test_lane, pAdjustMethod="BH")
test_lane_results <- test_lane_results[order(-test_lane_results$padj),]
head(test_lane_results)
summary(test_lane_results)
hist(na.omit(test_lane_results$pvalue))
plotDispEsts(test_lane, xlab="Mean of Normalized Counts", ylab="Dispersion", main="Mean Dispersion")
plotMA(test_lane, ylim=c(-3,3), main= "Diffential Expression by Lane")
summary(test_lane_results)


pca_analysis <- rlog(test_lane, blind=TRUE) 

#Testing for grouping by lane, age, or treatment
plotPCA(pca_analysis, intgroup=c("lane")) 
plotPCA(pca_analysis, intgroup=c("age"))
plotPCA(pca_analysis, intgroup=c("treatment"))


#Comparing each group with the control treatment 4 

#LPS_LPS vs vec_vec mice
DESeq_data_1_4 <- DESeqDataSetFromMatrix(tot_count_matrix, experimental_design, design = formula(~treatment + age + treatment:age))
DESeq_data_1_4$treatment <- relevel(DESeq_data_1_4$treatment, ref= "4")
DESeq_data_1_4 <- DESeq(DESeq_data_1_4)
resultsNames(DESeq_data_1_4)
treatment_results_1_vs_4 <- results(DESeq_data_1_4, contrast=list("treatment_1_vs_4", "age_Y_vs_O"), pAdjustMethod="BH" , alpha= 0.05)
summary(treatment_results_1_vs_4)
treatment_results_1_vs_4_sorted <- treatment_results_1_vs_4[order(treatment_results_1_vs_4$padj),]

write.csv(treatment_results_1_vs_4_sorted, file="treatment_diff_1_vs_4")
plotDispEsts(DESeq_data_1_4, xlab="Mean of Normalized Counts", ylab="Dispersion", main="Mean Dispersion")
plotMA(treatment_results_1_vs_4, ylim=c(-10,10), main="LPS LPS vs vec vec")

#LPS_vec vs vec_vec mice
DESeq_data_2_4 <- DESeqDataSetFromMatrix(tot_count_matrix, experimental_design, design = formula(~treatment + age + treatment:age))
DESeq_data_2_4$treatment <- relevel(DESeq_data_2_4$treatment, ref= "4")
DESeq_data_2_4 <- DESeq(DESeq_data_2_4)
resultsNames(DESeq_data_2_4)
treatment_results_2_vs_4 <- results(DESeq_data_2_4, contrast=list("treatment_2_vs_4", "age_Y_vs_O"), pAdjustMethod="BH" , alpha= 0.05)
summary(treatment_results_2_vs_4)
treatment_results_2_vs_4_sorted <- treatment_results_2_vs_4[order(treatment_results_2_vs_4$padj),]

write.csv(treatment_results_2_vs_4_sorted, file="treatment_diff_2_vs_4")
plotDispEsts(DESeq_data_2_4, xlab="Mean of Normalized Counts", ylab="Dispersion", main="Mean Dispersion")
plotMA(treatment_results_2_vs_4, ylim=c(-10,10), main="LPS vec vs vec vec")


#vec_LPS vs vec_vec mice
DESeq_data_3_4 <- DESeqDataSetFromMatrix(tot_count_matrix, experimental_design, design = formula(~treatment + age + treatment:age))
DESeq_data_3_4$treatment <- relevel(DESeq_data_3_4$treatment, ref= "4")
DESeq_data_3_4 <- DESeq(DESeq_data_3_4)
resultsNames(DESeq_data_3_4)
treatment_results_3_vs_4 <- results(DESeq_data_3_4, contrast=list("treatment_3_vs_4", "age_Y_vs_O"), pAdjustMethod="BH" , alpha= 0.05)
summary(treatment_results_3_vs_4)
treatment_results_3_vs_4_sorted <- treatment_results_3_vs_4[order(treatment_results_3_vs_4$padj),]

write.csv(treatment_results_3_vs_4_sorted, file="treatment_diff_3_vs_4")
plotDispEsts(DESeq_data_3_4, xlab="Mean of Normalized Counts", ylab="Dispersion", main="Mean Dispersion")
plotMA(treatment_results_3_vs_4, ylim=c(-10,10), main="vec LPS vs vec vec")


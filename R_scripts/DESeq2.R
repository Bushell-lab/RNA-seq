#load libraries
library(DESeq2)
library(tidyverse)
library(vsn)

Totals_filenames <- c('')

parent_dir <- ""

read_csv("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/FASTAs/human/GENCODE/v38/transcript_info/gencode.v38.pc_transcripts_gene_IDs.csv", col_names = c("transcript", "gene", "gene_sym"))  %>%
  group_by(gene) %>%
  sample_n(size = 1) -> gene_IDs

#read in fcounts----
data_list <- list()
for (sample in Totals_filenames) {
  read_tsv(file = file.path(parent_dir, "fCounts", paste0(sample, "_fCounts_output_matrix.txt")), skip = 2, col_names = c("gene", "counts")) %>%
    mutate(sample = rep(sample)) -> data_list[[sample]]
}

do.call("rbind", data_list) %>%
  spread(key = sample, value = counts) %>%
  filter(gene %in% gene_IDs$gene) %>%
  column_to_rownames("gene") -> counts_data
summary(counts_data)

#create a data frame with the condition/replicate information----
#you need to make sure this data frame is correct for your samples, the below creates one for a n=3 with EFT226 treatment.
sample_info <- data.frame(row.names = colnames(counts_data),
                          condition = factor(c(rep("CNOT", 3), rep("H3", 3), rep("Input", 3))),
                          replicate = factor(rep(1:3, 3)))

#print the data frame to visually check it has been made as expected
sample_info

#make a DESeq data set from imported data----
ddsTxi <- DESeqDataSetFromMatrix(counts_data,
                                   colData = sample_info,
                                   design = ~ condition + replicate)

#pre-filter to remove genes with less than an average of 10 counts across all samples----
keep <- rowMeans(counts(ddsTxi)) >= 10
table(keep)
ddsTxi <- ddsTxi[keep,]

#make sure levels are set appropriately so that Ctrl is "untreated"
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "Input")

#run DESeq on DESeq data set----
dds <- DESeq(ddsTxi)

#extract results for each comparison----
CNOT1_res <- results(dds, contrast=c("condition", "CNOT", "Input"))
H3_res <- results(dds, contrast=c("condition", "H3", "Input"))

#summarise results----
summary(CNOT1_res)
summary(H3_res)

#write summary to a text file
sink(file = file.path(parent_dir, "Analysis/DESeq2_output", "CNOT1_DEseq2_summary.txt"))
summary(CNOT1_res)
sink()

sink(file = file.path(parent_dir, "Analysis/DESeq2_output", "H3_DEseq2_summary.txt"))
summary(H3_res)
sink()

#apply LFC shrinkage for each comparison----
CNOT1_lfc_shrink <- lfcShrink(dds, coef="condition_CNOT_vs_Input", type="apeglm")
H3_lfc_shrink <- lfcShrink(dds, coef="condition_H3_vs_Input", type="apeglm")

#write reslts to csv----
as.data.frame(CNOT1_lfc_shrink[order(CNOT1_lfc_shrink$padj),]) %>%
  rownames_to_column("gene") %>%
  inner_join(gene_IDs, by = "gene") -> CNOT1_DEseq2_output
write_csv(CNOT1_DEseq2_output, file = file.path(parent_dir, "Analysis/DESeq2_output/CNOT1_DEseq2_apeglm_LFC_shrinkage.csv"))

as.data.frame(H3_lfc_shrink[order(H3_lfc_shrink$padj),]) %>%
  rownames_to_column("gene") %>%
  inner_join(gene_IDs, by = "gene") -> H3_DEseq2_output
write_csv(H3_DEseq2_output, file = file.path(parent_dir, "Analysis/DESeq2_output/H3_DEseq2_apeglm_LFC_shrinkage.csv"))

#extract normalised counts and plot SD vs mean----
ntd <- normTransform(dds) #this gives log2(n + 1)
vsd <- vst(dds, blind=FALSE) #Variance stabilizing transformation
rld <- rlog(dds, blind=FALSE) #Regularized log transformation

meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

#write out normalised counts data----
#Variance stabilizing transformation looks preferable for this data. Check for your own data and select the appropriate one
#The aim is for the range of standard deviations to be similar across the range of abundances, i.e. for the red line to be flat
as.data.frame(assay(vsd)) %>%
  rownames_to_column("gene") %>%
  inner_join(gene_IDs, by = "gene") -> normalised_counts
write_csv(normalised_counts, file = file.path(parent_dir, "Analysis/DESeq2_output/normalised_counts.csv"))

#plot PCA----
pcaData <- plotPCA(vsd, intgroup=c("condition", "replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png(filename = file.path(parent_dir, "plots/PCAs", "PCA.png"), width = 400, height = 350)
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
  geom_point(size=3) +
  #geom_text(aes(label=replicate), colour = 'black',size = 6, nudge_x = 2, vjust=1)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
dev.off()

#apply batch correct and re-plot heatmap and PCA----
mat <- assay(vsd)
mat <- limma::removeBatchEffect(mat, vsd$replicate)
assay(vsd) <- mat

#PCA
pcaData <- plotPCA(vsd, intgroup=c("condition", "replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png(filename = file.path(parent_dir, "plots/PCAs", "batch_corrected_PCA.png"), width = 400, height = 350)
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
  ggtitle("batch corrected")
dev.off()



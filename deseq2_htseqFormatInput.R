library(ggplot2)
library(factoextra)
library(scales)
library(writexl)
library("DESeq2")




#############
# Variables #
#############
plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots/wt6h_vs_ko6h/"
contrast=c("condition","KO_6", "WT_6")
selected_cond=c("KO_6", "WT_6")

design_path = "/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/experiment_design.txt"
count_path <- "/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/counts/"





################
# Data loading #
################
sampleFiles <- list.files(count_path)
sampleNames <- sub(".counts.qualityUniqueStranded.txt","",sampleFiles)



# Create a dictionary to transform the ENSMBL ID to Gene SYMBOL
genedictionary <-read.delim(paste0(count_path, sampleFiles[1]), stringsAsFactors = F, row.names = 1, header = F)[,"V2", drop=FALSE]



# Experiment design loading 
sampleConditionDic = read.table(design_path, header = FALSE, row.names = 1, stringsAsFactors = FALSE)
colnames(sampleConditionDic) = c("condition", "tiempo")
sampleConditionDic$pasted = paste(sampleConditionDic$condition, sampleConditionDic$tiempo, sep = "_")



# DESEq2 input table
sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleConditionDic[sampleNames,"condition"],
                          tiempo = sampleConditionDic[sampleNames,"tiempo"],
                          pasted = sampleConditionDic[sampleNames,"pasted"])

sampleTable$condition <- factor(sampleTable$condition)
sampleTable$tiempo <- factor(sampleTable$tiempo)
sampleTable$pasted <- factor(sampleTable$pasted)

sampleTable = sampleTable[sampleTable$pasted %in% selected_cond,]
sampleTable = sampleTable[,c("sampleName", "fileName", "pasted")]
colnames(sampleTable) = c("sampleName", "fileName", "condition")





##########
# DESeq2 #
##########

# DESeq2 data loading
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = count_path,
                                       design= ~ condition)



# Data filtering
keep <- rowSums(counts(ddsHTSeq)) >= 50
ddsHTSeq <- ddsHTSeq[keep,]

num_samples <- ncol(counts(ddsHTSeq))
sample_threshold <- round(num_samples*2/3)
keep <- apply(counts(ddsHTSeq), 1, function(x) as.integer(table(x)["0"]) < sample_threshold)
keep[is.na(keep)] <- TRUE
ddsHTSeq <- ddsHTSeq[keep,]



# DESeq2 Differential expression analysis
dds <- DESeq(ddsHTSeq)
resultsNames(dds)
res <- results(dds, contrast=contrast)
nrow(res)
normalizematrix = counts(dds, normalized=TRUE)



# Differencial expression analysis 
resdf = as.data.frame(res)
resdf = resdf[!is.na(resdf$baseMean),]
resdf$geneid = rownames(resdf) 
resdf$genename = genedictionary[rownames(resdf),"V2"]
write.table(resdf, paste0(plot_dir,"DESeq_results.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write_xlsx(resdf, paste0(plot_dir,"DESeq_results.xlsx"), col_names = TRUE)



# Genes which pass an adjusted p value threshold
res005 = subset(res, padj < 0.05)
res005df = data.frame(res005)
res005df$geneid = rownames(res005df) 
res005df$genename = genedictionary[rownames(res005df),"V2"]
write.table(res005df, paste0(plot_dir,"DESeq_results_0.05.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write_xlsx(res005df, paste0(plot_dir,"DESeq_results_0.05.xlsx"), col_names = T)
nrow(res005)





#########
# Plots #
#########
# Normalized count Distribution
boxplot(log2(normalizematrix+0.0001))
# Raw count Distribution
boxplot(log2(counts(dds, normalized=FALSE)+0.0001)) 
# Histogram of p-values from the call to nbinomTest.
hist(res$pvalue, breaks=100, col="skyblue", border="slateblue", main="")
# Histogram of p-adj from the call to nbinomTest.
hist(res$padj, breaks=100, col="darkorange1", border="darkorange4", main="")
# plotMA
jpeg(paste0(plot_dir,"MA_plot.jpeg"), width = 700, height = 500, quality = 100)
plotMA(res, main="DESeq2", ylim=c(-8,8))
dev.off()





# Normalized count Distribution
normalized_df = data.frame(sample = rep(colnames(normalizematrix), each=nrow(normalizematrix)),
                           values = as.vector(normalizematrix))
plot_norm = ggplot(normalized_df, aes(x = sample, y = values+0.0001)) +
  geom_boxplot() + 
  scale_y_continuous(trans='log2',breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 45)) +
  ylab("Counts") + xlab("Sample")
plot_norm
ggsave(path = plot_dir, filename = "Normaliced_counts_boxplot.jpeg", device = "jpeg", 
       plot = plot_norm, width = 7, height = 5, dpi = 500)


# Raw count Distribution
raw_df = data.frame(sample = rep(colnames(counts(dds, normalized=FALSE)), each=nrow(counts(dds, normalized=FALSE))),
                    values = as.vector(counts(dds, normalized=FALSE)))
plot_raw = ggplot(raw_df, aes(x = sample, y = values+0.0001)) +
  geom_boxplot() + 
  scale_y_continuous(trans='log2',breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 45)) +
  ylab("Counts") + xlab("Sample")
plot_raw
ggsave(path = plot_dir, filename = "Raw_counts_boxplot.jpeg", device = "jpeg", 
       plot = plot_raw, width = 7, height = 5, dpi = 500)





#######
# PCA #
#######
# rlogtranformation for PCA analysis
rld <- rlog(dds)
pca_plot = plotPCA(rld) + geom_text(aes(label = name), nudge_y = 0.5)
pca_plot

plotPCA(rld, returnData = TRUE)
write.table(plotPCA(rld, returnData = TRUE), paste0(plot_dir,"PCA_values.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")
write_xlsx(plotPCA(rld, returnData = TRUE), paste0(plot_dir,"PCA_values.xlsx"), col_names = T)

ggsave(path = plot_dir, filename = "PCA_DESeq2.jpeg", device = "jpeg", 
       plot = pca_plot, width = 7, height = 5, dpi = 500)



# library(pca3d)
# dds_matrix <- estimateSizeFactors(ddsHTSeq)
# 
# pca <- prcomp(t(counts(dds_matrix, normalized=TRUE)),scale.= TRUE)
# 
# pca3d(pca, components = 1:3, col = as.factor(c(rep("red",4), rep("blue",4))), title = NULL, new = TRUE, radius=1.4,
#       axes.color = "black", bg = "white", group = as.factor(c(rep("ApoE",4),rep("WT",4))), shape="sphere",
#       show.shadows = FALSE,show.plane = FALSE,
#       show.ellipses = TRUE, ellipse.ci = 0.65)





###############
# Marix count #
###############
dds_matrix <- estimateSizeFactors(ddsHTSeq)


conteoraw = data.frame(counts(dds_matrix))
conteoraw$geneid = rownames(conteoraw) 
conteoraw$genename = genedictionary[rownames(conteoraw),"V2"]
conteoraw = conteoraw[conteoraw$geneid %in% res005df$geneid, ]
write.table(conteoraw, paste0(plot_dir,"raw_counts_0.05.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write_xlsx(conteoraw, paste0(plot_dir,"raw_counts_0.05.xlsx"),col_names = T)

conteonorm = data.frame(counts(dds_matrix, normalized=TRUE),col_names = T)
conteonorm$geneid = rownames(conteonorm) 
conteonorm$genename = genedictionary[rownames(conteonorm),"V2"]
conteonorm = conteonorm[conteonorm$geneid %in% res005df$geneid, ]
write.table(conteonorm, paste0(plot_dir,"norm_counts_0.05.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write_xlsx(conteonorm, paste0(plot_dir,"norm_counts_0.05.xlsx"),col_names = T)


conteoraw = data.frame(counts(dds_matrix))
conteoraw$geneid = rownames(conteoraw) 
conteoraw$genename = genedictionary[rownames(conteoraw),"V2"]
write.table(conteoraw, paste0(plot_dir,"raw_counts_all.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write_xlsx(conteoraw, paste0(plot_dir,"raw_counts_all.xlsx"),col_names = T)

conteonorm = data.frame(counts(dds_matrix, normalized=TRUE))
conteonorm$geneid = rownames(conteonorm) 
conteonorm$genename = genedictionary[rownames(conteonorm),"V2"]
write.table(conteonorm, paste0(plot_dir,"norm_counts_all.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write_xlsx(conteonorm, paste0(plot_dir,"norm_counts_all.xlsx"),col_names = T)







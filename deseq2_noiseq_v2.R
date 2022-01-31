library(ggplot2)
library(factoextra)
library(scales)
library(writexl)



plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots/wt0h_vs_ko0h/"
contrast=c("condition","KO_0", "WT_0")
selected_cond=c("KO_0", "WT_0")

plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots/wt6h_vs_ko6h/"
contrast=c("condition","KO_6", "WT_6")
selected_cond=c("KO_6", "WT_6")

plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots/wt18h_vs_ko18h/"
contrast=c("condition","KO_18", "WT_18")
selected_cond=c("KO_18", "WT_18")



plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots/wt0h_vs_wt6h/"
contrast=c("condition","WT_6", "WT_0")
selected_cond=c("WT_6", "WT_0")

plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots/wt0h_vs_wt18h/"
contrast=c("condition","WT_18", "WT_0")
selected_cond=c("WT_18", "WT_0")



plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots/ko0h_vs_ko6h/"
contrast=c("condition","KO_6", "KO_0")
selected_cond=c("KO_6", "KO_0")

plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots/ko0h_vs_ko18h/"
contrast=c("condition","KO_18", "KO_0")
selected_cond=c("KO_18", "KO_0")



#plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots/temporal/"



################
# Data loading #
################


directory <- "/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/counts/"
sampleFiles <- list.files(directory)

genedictionary <-read.delim(paste0(directory, sampleFiles[1]), stringsAsFactors = F, row.names = 1, header = F)

sampleNames <- sub(".counts.qualityUniqueStranded.txt","",sampleFiles)


sampleConditionDic = read.table("/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/experiment_design.txt",header = FALSE, row.names = 1, stringsAsFactors = FALSE)
colnames(sampleConditionDic) = c("condition", "tiempo")
sampleConditionDic$pasted = paste(sampleConditionDic$condition, sampleConditionDic$tiempo, sep = "_")


sampleCondition <- sampleConditionDic[sampleNames,"condition"]
sampleTiempo <- sampleConditionDic[sampleNames,"tiempo"]
samplePasted <- sampleConditionDic[sampleNames,"pasted"]

sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleCondition,
                          tiempo = sampleTiempo,
                          pasted = samplePasted)

sampleTable$condition <- factor(sampleTable$condition)
sampleTable$tiempo <- factor(sampleTable$tiempo)
sampleTable$pasted <- factor(sampleTable$pasted)

sampleTable = sampleTable[sampleTable$pasted %in% selected_cond,]
sampleTable = sampleTable[,c("sampleName", "fileName", "pasted")]
colnames(sampleTable) = c("sampleName", "fileName", "condition")

##########
# DESeq2 #
##########

library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
ddsHTSeq




keep <- rowSums(counts(ddsHTSeq)) >= 50
ddsHTSeq <- ddsHTSeq[keep,]

num_samples <- ncol(counts(ddsHTSeq))
sample_threshold <- round(num_samples*2/3)
keep <- apply(counts(ddsHTSeq), 1, function(x) as.integer(table(x)["0"]) < sample_threshold)
keep[is.na(keep)] <- TRUE
ddsHTSeq <- ddsHTSeq[keep,]


dds <- DESeq(ddsHTSeq)
resultsNames(dds)
res <- results(dds, contrast=contrast)
res
nrow(res)

normalizematrix = counts(dds, normalized=TRUE)




# Distribution of count
boxplot(log2(normalizematrix+0.0001))
boxplot(log2(counts(dds, normalized=FALSE)+0.0001)) 


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







resdf = as.data.frame(res)
resdf = resdf[!is.na(resdf$baseMean),]
resdf$geneid = rownames(resdf) 
resdf$genename = genedictionary[rownames(resdf),"V2"]
# resdf$gen = mycountstmp[row.names(resdf),1]
# outdf = resdf[,c("gen", "log2FoldChange")]
write.table(resdf, paste0(plot_dir,"DESeq_results.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write_xlsx(resdf, paste0(plot_dir,"DESeq_results.xlsx"), col_names = TRUE)

##Histogram of p-values from the call to nbinomTest.
hist(res$pvalue, breaks=100, col="skyblue", border="slateblue", main="")

##Histogram of p-adj from the call to nbinomTest.
hist(res$padj, breaks=100, col="darkorange1", border="darkorange4", main="")


# plotMA
plotMA(res, main="DESeq2", ylim=c(-8,8))


# Genes which pass an adjusted p value threshold
res005 = subset(res, padj < 0.05)
res005df = data.frame(res005)
res005df$geneid = rownames(res005df) 
res005df$genename = genedictionary[rownames(res005df),"V2"]
write.table(res005df, paste0(plot_dir,"DESeq_results_0.05.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write_xlsx(res005df, paste0(plot_dir,"DESeq_results_0.05.xlsx"), col_names = T)
nrow(res005)


# rlogtranformation for PCA analysis
rld <- rlog(dds)
#head(assay(rld), 100)
plotPCA(rld) + geom_text(aes(label = name), nudge_y = 0.5)

plotPCA(rld, returnData = TRUE)
write.table(plotPCA(rld, returnData = TRUE), paste0(plot_dir,"PCA_values.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")
write_xlsx(plotPCA(rld, returnData = TRUE), paste0(plot_dir,"PCA_values.xlsx"), col_names = T)

ggsave(path = plot_dir, filename = "PCA_DESeq2.jpeg", device = "jpeg", 
       plot = plotPCA(rld) + geom_text(aes(label = name), nudge_y = 0.5) , width = 7, height = 5, dpi = 500)



library(pca3d)
dds_matrix <- estimateSizeFactors(ddsHTSeq)

pca <- prcomp(t(counts(dds_matrix, normalized=TRUE)),scale.= TRUE)

pca3d(pca, components = 1:3, col = as.factor(c(rep("red",4), rep("blue",4))), title = NULL, new = TRUE, radius=1.4,
      axes.color = "black", bg = "white", group = as.factor(c(rep("ApoE",4),rep("WT",4))), shape="sphere",
      show.shadows = FALSE,show.plane = FALSE,
      show.ellipses = TRUE, ellipse.ci = 0.65)






# matriz de conteo
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












##########
# NOISeq #
##########

# Data loading 
mycounts = NA
sampleNames = NULL
for (file in sampleFiles){
  fullname = paste(directory, file, sep = "")
  mycountstmp = read.table(fullname, header = FALSE, row.names = 1, fill = TRUE)
  sampleNames = c(sampleNames, sub(".counts.qualityUniqueStranded.txt", "", file))
  mycounts = cbind(mycounts, mycountstmp["V3"])
}
colnames(mycounts) = c("NA", sampleNames)
mycountstmp = mycounts[,sampleNames]
mycounts = na.omit(mycountstmp)
mycounts_copy = mycounts



# Sample filtering
plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots_NOIseq/wt0h_vs_ko0h/"
selected_cond=c("KO_0", "WT_0")

plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots_NOIseq/wt6h_vs_ko6h/"
selected_cond=c("KO_6", "WT_6")

plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots_NOIseq/wt18h_vs_ko18h/"
selected_cond=c("KO_18", "WT_18")



plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots_NOIseq/wt0h_vs_wt6h/"
selected_cond=c("WT_6", "WT_0")

plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots_NOIseq/wt0h_vs_wt18h/"
selected_cond=c("WT_18", "WT_0")



plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots_NOIseq/ko0h_vs_ko6h/"
selected_cond=c("KO_6", "KO_0")

plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots_NOIseq/ko0h_vs_ko18h/"
selected_cond=c("KO_18", "KO_0")


plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots_NOIseq//"
selected_cond=c("WT_6", "WT_0", "WT_18", "KO_18", "KO_0", "KO_6")




mycounts = mycounts_copy
sampleConditionDic = read.table("/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/experiment_design.txt",header = FALSE, row.names = 1, stringsAsFactors = FALSE)
colnames(sampleConditionDic) = c("condition", "tiempo")
sampleConditionDic$pasted = paste(sampleConditionDic$condition, sampleConditionDic$tiempo, sep = "_")
selected_samples = rownames(sampleConditionDic[sampleConditionDic$pasted %in% selected_cond, ])

mycounts = mycounts[,selected_samples]

mycounts = mycounts[rowSums(mycounts) >= 50, ]

num_samples <- ncol(mycounts)
sample_threshold <- round(num_samples*2/3)
keep <- apply(mycounts, 1, function(x) as.integer(table(x)["0"]) < sample_threshold)
keep[is.na(keep)] <- TRUE
mycounts <- mycounts[keep,]

myfactors = data.frame(condition = sampleConditionDic[selected_samples,"pasted"])


library(NOISeq)
mydata <- readData(data = mycounts, factors = myfactors)

# Quality control plots
mycountsbio = dat(mydata, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "boxplot")
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot")

mysaturation = dat(mydata, k = 0, ndepth = 7, type = "saturation")
explo.plot(mysaturation, toplot = 1, samples = 1:8)


# Data normalization
myRPKM = rpkm(assayData(mydata)$exprs, k = 0, lc = 1)
myUQUA = uqua(assayData(mydata)$exprs, lc = 0.5, k = 0)
myTMM = tmm(assayData(mydata)$exprs, long = 1000, lc = 0)
head(myRPKM[, 1:4])
head(myUQUA[, 1:4])
head(myTMM[, 1:4])

#sampleCondition <- sampleConditionDic[selected_samples,"pasted"]

sampleConditionDic_filtered = sampleConditionDic
sampleConditionDic_filtered$sampleName = row.names(sampleConditionDic_filtered)
#sampleConditionDic_filtered = sampleConditionDic_filtered[!sampleConditionDic_filtered$sampleName %in% c("FJD092", "FJD071", "CN023"),]
#sampleConditionDic_filtered = as.factor(sampleConditionDic_filtered[colnames(myTMM),"condition"])

mydata.norm = readData(data = myTMM, factors = sampleConditionDic_filtered[colnames(myTMM),])
mycountsbio = dat(mydata.norm, factor = NULL, type = "countsbio",  norm = TRUE, logtransf = FALSE)
explo.plot(mycountsbio, toplot = "global", plottype = "boxplot", samples = 1:ncol(myTMM))
explo.plot(mycountsbio, toplot = "global", samples = 1:ncol(myTMM), plottype = "barplot")
res.pca <- prcomp(t(log2(exprs(mydata.norm)+0.0001)), scale. = TRUE)

plot(res.pca$rotation[,1], res.pca$rotation[,2])



fviz_eig(res.pca)

pca_noiseq = fviz_pca_ind(res.pca, axes = c(1, 2), 
             col.ind = sampleConditionDic[colnames(exprs(mydata.norm)),"pasted"], # Color by the quality of representation
             repel = TRUE # Avoid text overlapping
             )
pca_noiseq
ggsave(path = plot_dir, filename = "PCA_NOISeq.jpeg", device = "jpeg", 
       plot = pca_noiseq, width = 7, height = 5, dpi = 500)

ggsave(path = plot_dir, filename = "PCA_NOISeq.pdf", device = "pdf", 
       plot = pca_noiseq, width = 7, height = 5, dpi = 500)


# Differential expression
myresults <- noiseq(mydata, factor = "condition", k = NULL, norm = "tmm", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "no")



mynoiseq.deg = degenes(myresults, q = 0.9, M = NULL)
mynoiseq.deg1 = degenes(myresults, q = 0.9, M = "up")
mynoiseq.deg2 = degenes(myresults, q = 0.9, M = "down")

DE.plot(myresults, q = 0.9, graphic = "expr", log.scale = TRUE)



mynoiseq.deg$geneid = rownames(mynoiseq.deg) 
mynoiseq.deg$genename = genedictionary[rownames(mynoiseq.deg),"V2"]
write.table(mynoiseq.deg, paste0(plot_dir,"NOIseq_results_0.9.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write_xlsx(mynoiseq.deg, paste0(plot_dir,"NOIseq_results_0.9.xlsx"), col_names = T)
dim(mynoiseq.deg)



# Matix count

raw_counts = data.frame(assayData(mydata)$exprs, stringsAsFactors = F)
raw_counts$geneid = rownames(raw_counts)
raw_counts$genename = genedictionary[rownames(raw_counts),"V2"]


normalized_counts = data.frame(myTMM, stringsAsFactors = F)
normalized_counts$geneid = rownames(normalized_counts)
normalized_counts$genename = genedictionary[rownames(normalized_counts),"V2"]




write.table(raw_counts, paste0(plot_dir,"raw_counts_all.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write_xlsx(raw_counts, paste0(plot_dir,"raw_counts_all.xlsx"),col_names = T)


write.table(normalized_counts, paste0(plot_dir,"norm_counts_all.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write_xlsx(normalized_counts, paste0(plot_dir,"norm_counts_all.xlsx"),col_names = T)


raw_counts80 = raw_counts[mynoiseq.deg$geneid,]
write.table(raw_counts80, paste0(plot_dir,"raw_counts_90.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write_xlsx(raw_counts80, paste0(plot_dir,"raw_counts_90.xlsx"), col_names = T)


normalized_counts80 = normalized_counts[mynoiseq.deg$geneid,]
write.table(normalized_counts80, paste0(plot_dir,"norm_counts_90.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write_xlsx(normalized_counts80, paste0(plot_dir,"norm_counts_90.xlsx"), col_names = T)









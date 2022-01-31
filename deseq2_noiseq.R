library(ggplot2)
library(factoextra)
library(scales)


plot_dir="/home/gonzalo/UAMssh/fjd/BioinfoUnit/NataliaMarta_Vito/Marta_reagrupacion/plots/wt_vs_ko/"
plot_dir="/home/gonzalo/UAMssh/fjd/BioinfoUnit/NataliaMarta_Vito/Marta_reagrupacion/plots/wt_vs_wt6h/"
plot_dir="/home/gonzalo/UAMssh/fjd/BioinfoUnit/NataliaMarta_Vito/Marta_reagrupacion/plots/wt_vs_wt18h/"
plot_dir="/home/gonzalo/UAMssh/fjd/BioinfoUnit/NataliaMarta_Vito/Marta_reagrupacion/plots/wt_vs_ko6h/"
plot_dir="/home/gonzalo/UAMssh/fjd/BioinfoUnit/NataliaMarta_Vito/Marta_reagrupacion/plots/wt_vs_ko18h/"
plot_dir="/home/gonzalo/UAMssh/fjd/BioinfoUnit/NataliaMarta_Vito/Marta_reagrupacion/plots/ko_vs_ko6h/"
plot_dir="/home/gonzalo/UAMssh/fjd/BioinfoUnit/NataliaMarta_Vito/Marta_reagrupacion/plots/ko_vs_ko18h/"
plot_dir="/home/gonzalo/UAMssh/fjd/BioinfoUnit/NataliaMarta_Vito/Marta_reagrupacion/plots/temporal/"


################
# Data loading #
################
directory <- "/home/gonzalo/UAMssh/fjd/BioinfoUnit/NataliaMarta_Vito/Marta_reagrupacion/counts/"
sampleFiles <- list.files(directory)


sampleNames <- sub("_.*counts.qualityUniqueStranded.txt","",sampleFiles)


sampleConditionDic = read.table("/home/gonzalo/UAMssh/fjd/BioinfoUnit/NataliaMarta_Vito/Marta_reagrupacion/experiment_design.tsv",header = FALSE, row.names = 1, stringsAsFactors = FALSE)
colnames(sampleConditionDic) = c("condition")


sampleCondition <- sampleConditionDic[sampleNames,"condition"]

sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleCondition)

sampleTable$condition <- factor(sampleTable$condition)





##########
# DESeq2 #
##########

library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
ddsHTSeq


dds <- DESeq(ddsHTSeq)
resultsNames(dds)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

res <- results(dds, contrast=c("condition","LATE_ONSET","EARLY_ONSET"))
res

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
  theme(axis.text.x = element_text(angle = 45)) +
  ylab("Counts") + xlab("Sample")
plot_norm
ggsave(path = plot_dir, filename = "Normaliced_counts_boxplot.jpeg", device = "jpeg", 
       plot = plot_norm, width = 10, height = 5, dpi = 500)



raw_df = data.frame(sample = rep(colnames(counts(dds, normalized=FALSE)), each=nrow(counts(dds, normalized=FALSE))),
                           values = as.vector(counts(dds, normalized=FALSE)))
plot_raw = ggplot(raw_df, aes(x = sample, y = values+0.0001)) +
  geom_boxplot() + 
  scale_y_continuous(trans='log2',breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  theme(axis.text.x = element_text(angle = 45)) +
  ylab("Counts") + xlab("Sample")
plot_raw
ggsave(path = plot_dir, filename = "Raw_counts_boxplot.jpeg", device = "jpeg", 
       plot = plot_raw, width = 10, height = 5, dpi = 500)






nrow(res)

resdf = as.data.frame(res)
resdf$miRNA = rownames(resdf) 
# resdf$gen = mycountstmp[row.names(resdf),1]
# outdf = resdf[,c("gen", "log2FoldChange")]
write.table(resdf, paste(plot_dir,"DESeq_results.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


##Histogram of p-values from the call to nbinomTest.
hist(res$pvalue, breaks=100, col="skyblue", border="slateblue", main="")

##Histogram of p-adj from the call to nbinomTest.
hist(res$padj, breaks=100, col="darkorange1", border="darkorange4", main="")


# plotMA
plotMA(res, main="DESeq2", ylim=c(-8,8))


# Genes which pass an adjusted p value threshold
res005 = subset(res, padj < 0.05)

# rlogtranformation for PCA analysis
rld <- rlog(dds)
head(assay(rld), 100)
plotPCA(rld) + geom_text(aes(label = name), nudge_y = 0.5)

plotPCA(rld, returnData = TRUE)

ggsave(path = plot_dir, filename = "PCA_DESeq2.jpeg", device = "jpeg", 
       plot = plotPCA(rld) + geom_text(aes(label = name), nudge_y = 0.5), width = 7, height = 5, dpi = 500)






##########
# NOISeq #
##########

# Data loading 
mycounts = NA
sampleNames = NULL
for (file in sampleFiles){
  fullname = paste(directory, file, sep = "")
  mycountstmp = read.table(fullname, header = FALSE, row.names = 1, fill = TRUE)
  sampleNames = c(sampleNames, sub("_.*", "", file))
  mycounts = cbind(mycounts, mycountstmp["V3"])
}
colnames(mycounts) = c("NA", sampleNames)
mycountstmp = mycounts[,sampleNames]
mycounts = na.omit(mycountstmp) 
mycounts = mycounts[rowSums(mycounts) >= 10, ]

myfactors = data.frame(condition = sampleCondition)

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


mydata.norm = readData(data = myTMM, factors = sampleConditionDic)
mycountsbio = dat(mydata.norm, factor = NULL, type = "countsbio",  norm = TRUE, logtransf = FALSE)
explo.plot(mycountsbio, toplot = "global", plottype = "boxplot", samples = 1:nrow(sampleConditionDic))
explo.plot(mycountsbio, toplot = "global", samples = 1:nrow(sampleConditionDic), plottype = "barplot")
res.pca <- prcomp(t(log2(exprs(mydata.norm)+0.0001)), scale. = TRUE)

plot(res.pca$rotation[,1], res.pca$rotation[,2])



fviz_eig(res.pca)

pca_noiseq = fviz_pca_ind(res.pca, axes = c(1, 2), 
             col.ind = sampleConditionDic[colnames(exprs(mydata.norm)),], # Color by the quality of representation
             repel = TRUE # Avoid text overlapping
             )
pca_noiseq
ggsave(path = plot_dir, filename = "PCA_NOISeq.jpeg", device = "jpeg", 
       plot = pca_noiseq, width = 7, height = 5, dpi = 500)




# Differential expression
myresults <- noiseq(mydata, factor = "condition", k = NULL, norm = "tmm", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "no")



mynoiseq.deg = degenes(myresults, q = 0.8, M = NULL)
mynoiseq.deg1 = degenes(myresults, q = 0.8, M = "up")
mynoiseq.deg2 = degenes(myresults, q = 0.8, M = "down")

DE.plot(myresults, q = 0.9, graphic = "expr", log.scale = TRUE)









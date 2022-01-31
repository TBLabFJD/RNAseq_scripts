library(ggplot2)
library(factoextra)
library(scales)
library(ggrepel)


plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots/wt_vs_ko/"
plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots/wt_vs_wt6h/"
plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots/wt_vs_wt18h/"
plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots/wt_vs_ko6h/"
plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots/wt_vs_ko18h/"
plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots/ko_vs_ko6h/"
plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots/ko_vs_ko18h/"
plot_dir="/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/plots/temporal/"


setwd(plot_dir)
################
# Data loading #
################

# Matrix count loading 
count_directory <- "/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/counts/"
sampleFiles <- list.files(count_directory)
mycounts = NA
sampleNames = NULL
for (file in sampleFiles){
  fullname = paste(count_directory, file, sep = "")
  mycountstmp = read.table(fullname, header = FALSE, row.names = 1, fill = TRUE)
  sampleNames = c(sampleNames, sub(".counts.qualityUniqueStranded.txt", "", file))
  mycounts = cbind(mycounts, mycountstmp["V3"])
}
colnames(mycounts) = c("NA", sampleNames)
mycountstmp = mycounts[,sampleNames]
mycounts = na.omit(mycountstmp) 



# Experiment design loading 
design_directory <- "/home/gonzalo/UAMssh/NOBACKUP/BioinfoUnit/LuisBlanco/RNAseq2/experiment_design.txt"
sampleConditionDic = read.table(design_directory,header = FALSE, row.names = 1, stringsAsFactors = FALSE)
colnames(sampleConditionDic) = c("condition", "tiempo")

sampleCondition <- sampleConditionDic[sampleNames,] # It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order.



# miRNA filtering
mycounts = mycounts[rowSums(mycounts) >= 50, ]

num_samples <- ncol(mycounts)
sample_threshold <- round(num_samples*2/3)
keep <- apply(mycounts, 1, function(x) as.integer(table(x)["0"]) < sample_threshold)
keep[is.na(keep)] <- TRUE
mycounts <- mycounts[keep,]




# Subset
mycounts_tmp = mycounts[,-which(names(mycounts) %in% c("41_4w"))]
sampleCondition_tmp = sampleCondition[-which(row.names(sampleCondition) %in% c("41_4w")),]



# # Data loading #
# 
# directory <- "/home/gonzalo/UAMssh/fjd/BioinfoUnit/miRNA_Celia_Perales/celia_experimento_3EV/counts/"
# sampleFiles <- list.files(directory)
# 
# 
# sampleNames <- sub("_R1_all_l17_25.counts.qualityUniqueStranded.txt","",sampleFiles)
# 
# 
# sampleConditionDic = read.table("/home/gonzalo/UAMssh/fjd/BioinfoUnit/miRNA_Celia_Perales/celia_experimento_3EV/experiment_design.txt",header = FALSE, row.names = 1, stringsAsFactors = FALSE)
# colnames(sampleConditionDic) = c("cirrosis", "mejoria", "tiempo")
# 
# 
# sampleCondition <- sampleConditionDic[sampleNames,c("cirrosis", "mejoria", "tiempo")]
# 
# sampleTable <- data.frame(sampleName = sampleNames,
#                           fileName = sampleFiles,
#                           sampleCondition)
# 
# sampleTable$cirrosis <- factor(sampleTable$cirrosis)
# sampleTable$mejoria <- factor(sampleTable$mejoria)
# #sampleTable$tiempo <- factor(sampleTable$tiempo, levels = c("Bas", "4w", "12w", "12wP", "48wP", "2aP"))
# sampleTable$tiempo <- time_dict[sampleTable$tiempo,"values"]
# 
# 
# 
# #sampleTable = sampleTable[!sampleTable$sampleName %in% c("41_4w", "155_2aP", "41_2aP"),]
# sampleTable = sampleTable[!sampleTable$sampleName %in% c("41_4w"),]





##########
# DESeq2 #
##########

####
# Data loading and filtering
####


library("DESeq2")
ddsHTSeq <- DESeqDataSetFromMatrix(countData = mycounts_tmp,
                                   colData = sampleCondition_tmp,
                                   design= ~ cirrosis + mejoria + tiempo)

# ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
#                                        directory = directory,
#                                        design= ~ cirrosis + mejoria + tiempo + cirrosis:tiempo + mejoria:tiempo)
ddsHTSeq


dds <- DESeq(ddsHTSeq)
resultsNames(dds)

keep <- rowSums(counts(dds)) >= 50
table(keep)
dds <- dds[keep,]

num_samples <- ncol(counts(dds))
sample_threshold <- round(num_samples*2/3)
keep <- apply(counts(dds), 1, function(x) as.integer(table(x)["0"]) < sample_threshold)
keep[is.na(keep)] <- TRUE
table(keep)
dds <- dds[keep,]



####
# Distribution of count
####
normalizematrix = counts(dds, normalized=TRUE)

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




####
# Differential expression
####
# Cirrótico vs no cirrótico
res <- results(dds, contrast=c("cirrosis","cirrotico","no_cirrotico"))

jpeg("LogFoldChange_cirrosis.jpeg", width = 700, height = 500, quality = 100)
plotMA(res, main="DESeq2", ylim=c(-8,8))
dev.off()

resdf = as.data.frame(res)
resdf$miRNA = rownames(resdf)
write.table(resdf, paste(plot_dir,"DESeq_results_cirroticos_vs_nocirroticos.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


# Mejora vs Empeora
res <- results(dds, contrast=c("mejoria","mejora","empeora"))

jpeg("LogFoldChange_mejoria.jpeg", width = 700, height = 500, quality = 100)
plotMA(res, main="DESeq2", ylim=c(-8,8))
dev.off()

resdf = as.data.frame(res)
resdf$miRNA = rownames(resdf)
write.table(resdf, paste(plot_dir,"DESeq_results_mejora_vs_empeora.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")




# resdf = as.data.frame(res)
# resdf$miRNA = rownames(resdf) 
# resdf$gen = mycountstmp[row.names(resdf),1]
# outdf = resdf[,c("gen", "log2FoldChange")]
# write.table(resdf, paste(plot_dir,"DESeq_results.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# # Genes which pass an adjusted p value threshold
# res005 = subset(res, padj < 0.05)
# res005 = as.data.frame(res005)
# nrow(res005)


# # plotMA
# jpeg("LogFoldChange.jpg", width = 750, height = 350)
# plotMA(res, main="DESeq2", ylim=c(-8,8))
# dev.off()



# rlogtranformation for PCA analysis
rld <- rlog(dds)
head(assay(rld), 10)
plotPCA(rld, intgroup = "cirrosis") + geom_text_repel(aes(label = name), size = 3)
plotPCA(rld, intgroup = "mejoria") + geom_text_repel(aes(label = name), size = 3)
plotPCA(rld, intgroup = "tiempo") + geom_text_repel(aes(label = name), size = 3)



ggsave(path = plot_dir, filename = "PCA_DESeq2_cirrosis.jpeg", device = "jpeg", 
       plot = plotPCA(rld, intgroup = "cirrosis") + geom_text_repel(aes(label = name), size = 3), width = 7, height = 5, dpi = 500)

ggsave(path = plot_dir, filename = "PCA_DESeq2_mejoria.jpeg", device = "jpeg", 
       plot = plotPCA(rld, intgroup = "mejoria") + geom_text_repel(aes(label = name), size = 3), width = 7, height = 5, dpi = 500)

ggsave(path = plot_dir, filename = "PCA_DESeq2_tiempo.jpeg", device = "jpeg", 
       plot = plotPCA(rld, intgroup = "tiempo") + geom_text_repel(aes(label = name), size = 3), width = 7, height = 5, dpi = 500)

#design= ~ cirrosis + mejoria + tiempo + cirrosis:tiempo + mejoria:tiempo)





# 
# 
# ###############
# # Time series #
# ###############
# 
# library("DESeq2")
# # Cirrosis
# ddsHTSeq <- DESeqDataSetFromMatrix(countData = mycounts,
#                                    colData = sampleCondition,
#                                    design = ~ cirrosis  + tiempo + cirrosis:tiempo)
# ddsHTSeq
# ddsTC <- DESeq(ddsHTSeq, test="LRT", reduced = ~ cirrosis + tiempo)
# resTC <- results(ddsTC)
# resTC$symbol <- mcols(ddsTC)$symbol
# head(resTC[order(resTC$padj),], 4)
# 
# 
# mirna = row.names(resTC[resTC$padj < 0.1,])
# 
# 
# fiss <- plotCounts(ddsTC,mirna[4], 
#                    intgroup = c("tiempo","cirrosis"), returnData = TRUE)
# fiss$tiempo <- as.numeric(as.character(fiss$tiempo))
# ggplot(fiss,
#        aes(x = tiempo, y = count, color = cirrosis, group = cirrosis)) + 
#   geom_point() + stat_summary(fun=mean, geom="line") +
#   scale_y_log10()
# 
# 
# 
# 
# # Mejoría
# ddsHTSeq <- DESeqDataSetFromMatrix(countData = mycounts,
#                                    colData = sampleCondition,
#                                    design = ~ mejoria + tiempo + mejoria:tiempo)
# ddsHTSeq
# ddsTC <- DESeq(ddsHTSeq, test="LRT", reduced = ~ mejoria + tiempo)
# resTC <- results(ddsTC)
# resTC$symbol <- mcols(ddsTC)$symbol
# head(resTC[order(resTC$padj),], 4)
# 
# 
# mirna = row.names(resTC[resTC$padj < 0.05,])
# 
# 
# fiss <- plotCounts(ddsTC,mirna[1], 
#                    intgroup = c("tiempo","cirrosis"), returnData = TRUE)
# fiss$tiempo <- as.numeric(as.character(fiss$tiempo))
# ggplot(fiss,
#        aes(x = tiempo, y = count, color = cirrosis, group = cirrosis)) + 
#   geom_point() + stat_summary(fun=mean, geom="line") +
#   scale_y_log10()












###############
# Time series #
###############
library(maSigPro)
library(NOISeq)


# Data mormalization
myTMM = tmm(mycounts, long = 1000, lc = 0)


####
# Distribution of count
####
# Raw


raw_df = data.frame(sample = rep(colnames(mycounts), each=nrow(mycounts)),
                    values = as.vector(as.matrix(mycounts)))
plot_raw = ggplot(raw_df, aes(x = sample, y = values+0.0001)) +
  geom_boxplot() + 
  scale_y_continuous(trans='log2',breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  theme(axis.text.x = element_text(angle = 90, hjust=0.95)) +
  ylab("Counts") + xlab("Sample")
plot_raw
ggsave(path = plot_dir, filename = "Raw_counts_boxplot_all_samples.jpeg", device = "jpeg", 
       plot = plot_raw, width = 10, height = 5, dpi = 500)


# Normalized
normalized_df = data.frame(sample = rep(colnames(myTMM), each=nrow(myTMM)),
                           values = as.vector(myTMM))
plot_norm = ggplot(normalized_df, aes(x = sample, y = values+0.0001)) +
  geom_boxplot() + 
  scale_y_continuous(trans='log2',breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  theme(axis.text.x = element_text(angle = 90, hjust=0.95)) +
  ylab("Counts") + xlab("Sample")
plot_norm
ggsave(path = plot_dir, filename = "Normaliced_counts_boxplot_all_samples.jpeg", device = "jpeg", 
       plot = plot_norm, width = 10, height = 5, dpi = 500)








# Design creation

myfactors = sampleCondition
myfactors$Time = myfactors$tiempo

myfactors$KO = myfactors$condition
myfactors$KO = gsub("KO", 1, myfactors$KO)
myfactors$KO = gsub("WT", 0, myfactors$KO)
myfactors$WT = myfactors$condition 
myfactors$WT = gsub("KO", 0, myfactors$WT)
myfactors$WT = gsub("WT", 1, myfactors$WT)



myfactors$Replicates = as.numeric(as.factor(paste0(myfactors$condition, myfactors$tiempo)))


#####
# Category selection
#####

myfactors_sub = myfactors[,c("Time","Replicates","KO", "WT")]

#



# Create a regression matrix for the full regression model:
for (i in 1:ncol(myfactors_sub)) myfactors_sub[,i] = as.integer(myfactors_sub[,i])
design <- make.design.matrix(myfactors_sub, degree = 2)

myTMM_short = myTMM[1:8000,]
myTMM_short = myTMM[9229:9229,,drop = F]

selectedcol = myTMM[,which(grepl("C", colnames(myTMM)))]
suma_col = rowSums(selectedcol)
zeroC = which(suma_col == 0)
zeroCdf = myTMM[zeroC,]
View(zeroCdf)

zerolist = apply(myTMM, 1, function(x) as.numeric(table(x)["0"]))



myTMM_short = myTMM[9229:9230,,drop = F]
myTMM_short[1,2] == 6.351263
myTMM_short[1,2] = as.numeric(myTMM_short[1,2])
myTMM_short[1,2] = 6.351263

myTMM_short = myTMM[12647:12650,]
myTMM_short[1,2] = 5.08101



myTMM_short = myTMM
myTMM_short[9229,2] <- 6.351263
myTMM_short[12647,2] <- 5.08101



# Finding significant genes
fit <- p.vector(myTMM_short, design, Q = 0.05, min.obs = 20)

fit$i
fit$alfa
# fit$SELEC

# Finding significant differences
tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)


# Obtaining lists of significant genes
sigs_all <- get.siggenes(tstep, rsq = 0.7, vars = "all")
sigs_groups <- get.siggenes(tstep, rsq = 0.7, vars = "groups")
sigs_each <- get.siggenes(tstep, rsq = 0.7, vars = "each")

sigs_all$summary
sigs_groups$summary
sigs_each$summary


#suma2Venn(sigs_all$summary)
suma2Venn(sigs_groups$summary)
#suma2Venn(sigs_each$summary[, c(2,3,4,5)])

setdiff(sigs_groups$summary$WTvsKO, sigs_groups$summary$KO)

mycluster = see.genes(sigs_groups$sig.genes$WTvsKO, k = 9)
mycluster_geneName = mycluster$cut 
cluster_list = split(names(mycluster_geneName), mycluster_geneName)


mycluster_geneName_df = data.frame(mycluster_geneName)
mycluster_geneName_df$gene = row.names(mycluster_geneName_df)
colnames(mycluster_geneName_df) = c("cluster", "gene")
write.table(mycluster_geneName_df, "genes_cluster_temporal.txt", col.names = T, row.names = F, quote = F, sep = "\t")



see.genes(sigs_groups$sig.genes$WTvsKO, show.fit = T, dis =design$dis,
          cluster.method="hclust" ,cluster.data = 1, k = 7)



see.genes(sigs$sig.genes$mejora, show.fit = T, dis =design$dis,
          cluster.method="hclust" ,cluster.data = 1, k = 9)

see.genes(sigs$sig.genes$cirrotico, show.fit = T, dis =design$dis,
          cluster.method="hclust" ,cluster.data = 1, k = 4)







##########
# NOISeq #
##########
library(NOISeq)


mydata <- readData(data = mycounts_tmp, factors = sampleCondition_tmp)

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

sampleCondition <- sampleConditionDic[sampleNames,"cirrosis"]
sampleCondition <- sampleConditionDic[sampleNames,"mejoria"]


mydata.norm = readData(data = myTMM, factors = sampleCondition_tmp[colnames(myTMM),])
mycountsbio = dat(mydata.norm, factor = NULL, type = "countsbio",  norm = TRUE, logtransf = FALSE)
explo.plot(mycountsbio, toplot = "global", plottype = "boxplot", samples = 1:nrow(sampleCondition_tmp))
explo.plot(mycountsbio, toplot = "global", samples = 1:nrow(sampleCondition_tmp), plottype = "barplot")
res.pca <- prcomp(t(log2(exprs(mydata.norm)+0.0001)), scale. = TRUE)

plot(res.pca$rotation[,1], res.pca$rotation[,2])



fviz_eig(res.pca)


pca_noiseq = fviz_pca_ind(res.pca, axes = c(1, 2), 
                        col.ind = sampleConditionDic[colnames(exprs(mydata.norm)),"cirrosis"], # Color by the quality of representation
                        repel = TRUE, # Avoid text overlapping
                        show.legend = FALSE) + labs(color = "cirrosis", shape = "cirrosis") 
pca_noiseq
ggsave(path = plot_dir, filename = "PCA_NOISeq_cirrosis.jpeg", device = "jpeg", 
       plot = pca_noiseq, width = 7, height = 5, dpi = 500)


pca_noiseq = fviz_pca_ind(res.pca, axes = c(1, 2), 
                          col.ind = sampleConditionDic[colnames(exprs(mydata.norm)),"mejoria"], # Color by the quality of representation
                          repel = TRUE, # Avoid text overlapping
                          show.legend = FALSE) + labs(color = "mejoria", shape = "mejoria") 
pca_noiseq
ggsave(path = plot_dir, filename = "PCA_NOISeq_mejoria.jpeg", device = "jpeg", 
       plot = pca_noiseq, width = 7, height = 5, dpi = 500)




# Differential expression
myresults <- noiseq(mydata, factor = "condition", k = NULL, norm = "tmm", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "no")



mynoiseq.deg = degenes(myresults, q = 0.8, M = NULL)
mynoiseq.deg1 = degenes(myresults, q = 0.8, M = "up")
mynoiseq.deg2 = degenes(myresults, q = 0.8, M = "down")

DE.plot(myresults, q = 0.9, graphic = "expr", log.scale = TRUE)










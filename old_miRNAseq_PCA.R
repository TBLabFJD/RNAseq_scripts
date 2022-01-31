
# Date: 16/03/2020

library(readxl)
library(NOISeq)
library(factoextra)
library(OmicsPLS)
library(mixOmics)


#########################################
# Analysis from filtered Cufflinks data #
#########################################



# READ QUANTIFICATION

counts <- data.frame(read_excel("/Users/luisfelipeguerrero/BioinfoUnit/miRNA_Celia_Perales/genes_datosfiltrados-miRNAs-191210.xlsx", 
                                       sheet = "ALL TOGHETER"), stringsAsFactors = F)

mirna.counts <- counts[grepl("MIR", x = counts$Gene.name),]
rownames(mirna.counts) <- mirna.counts$Gene.name
mirna.counts <- mirna.counts[,3:(ncol(mirna.counts)-1)]

myfactors = data.frame(
  row.names=colnames(mirna.counts), 
  library=sapply(colnames(mirna.counts), function(x) strsplit(x,split = ".", fixed = T)[[1]][1])
)


# PCA

## Library load
library(openxlsx)
library(NOISeq)
library(factoextra)
library(OmicsPLS)
library(mixOmics)
library(DESeq2)


#mirna.counts.filtered = filtered.data(mymatrix[,cols], factor=as.factor(myfactors$condition), norm = FALSE,  method = 1, cv.cutoff = 100, cpm = 1, p.adj = "fdr")


seqDepth=t(as.data.frame(colSums(mirna.counts)))
barplot(seqDepth[1,], las=2)


## create dataset 

mydata = readData(data = mirna.counts, factors = myfactors)


########
## QC ##
########

## Count distribution per biotype

mycountsbio = dat(mydata, factor = NULL, type = "countsbio",  norm = TRUE, logtransf = FALSE)
explo.plot(mycountsbio, toplot = "global", plottype = "boxplot", samples = 1:nrow(myfactors))
explo.plot(mycountsbio, toplot = "global", samples = 1:nrow(myfactors), plottype = "barplot")


# RNA composition

mycomp = dat(mydata, norm = TRUE, logtransf = FALSE, type = "cd")
explo.plot(mycomp, samples = 1:nrow(myfactors))



# PCA

myPCA = dat(mydata, type = "PCA", norm = TRUE, logtransf = FALSE)
explo.plot(myPCA, factor = "library")

res.pca <- prcomp(t(log2(exprs(mydata)+0.0001)), scale. = FALSE)

plot(res.pca$rotation[,1], res.pca$rotation[,2])

fviz_eig(res.pca)

fviz_pca_ind(res.pca, axes = c(1, 2), 
             col.ind = myfactors$library, # Color by the quality of representation
             repel = TRUE # Avoid text overlapping
             
)

fviz_pca_ind(res.pca, axes = c(1, 3), 
             col.ind = myfactors$library, # Color by the quality of representation
             repel = TRUE # Avoid text overlapping
             
)


fviz_pca_biplot(res.pca, label="var",
                select.ind = list(contrib = 2))

####################################
# Analysis from raw Cufflinks data #
####################################



# READ QUANTIFICATION









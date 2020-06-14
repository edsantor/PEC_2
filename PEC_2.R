## ----setup, include=FALSE--------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----echo=TRUE, message=FALSE, warning=FALSE, results='hide'---------------------
library(readr)
counts = read_delim("datos/counts.csv", 
    ";", escape_double = FALSE, trim_ws = TRUE)


## ----echo=TRUE, message=FALSE, warning=FALSE, results='hide'---------------------
counts$X1 = gsub("\\..*", "", counts$X1)


## ----echo=TRUE, message=FALSE, warning=FALSE, results='hide'---------------------
targets = read_csv("datos/targets.csv")


## ----echo=TRUE, message=FALSE, warning=FALSE, results='hide'---------------------
# Establece una semilla para asegurar la repetitividad aleatoria.
set.seed(73026030)
# Crea un data frame vacío que contendrá los nombres de las 30 muestras.
df = data.frame()
# Almacena en groups los diferentes grupos de la columna Group del
# data frame targets.
groups = unique(targets[,"Group"])
# Bucle for que se ejecuta tantas veces como grupos haya en groups.
for (i in 1:nrow(groups)) {
  # Crea un data frame con únicamente las filas correspondientes al 
  # grupo i.
  targets_group = subset(targets, Group == as.character(groups[i,1]))
  # Cuenta las filas del data frame targets_groups, es decir, el número 
  # de muestras del grupo i.
  rows = nrow(targets_group)
  # Selecciona aleatoriamente 10 índices de fila sin repetición entre 
  # todas las filas de targets_groups del grupo i.
  indexes = sample(1:rows,10,replace = F)
  # Crea el data frame targets_samples que contiene los nombres de las 
  # muestras correspondientes a los 10 índices aleatorios seleccionados.
  targets_samples = targets_group[indexes,"Sample_Name"]
  # En cada iteración del bucle for, almacena en df los nombres de las 10
  # muestras de cada grupo que se han seleccionado aleatoriamnte.
  df = rbind(df,targets_samples)
}
# Crea el vector de caracteres df a partir del data frame df con los 
# nombres de las 30 muestras seleccionadas.
df = df$Sample_Name
# Crea la matriz de conteo de las 30 muestras, indenxando las columnas de 
# counts con el vector df. 
counts_set = counts[,df]
# Ordena alfanuméricamente las columnas de la matriz de conteo counts_set.
counts_set = counts_set[,order(colnames(counts_set))]
# Establece que los nombres de fila de counts_set sean los ID de Ensembl 
# de counts.
row.names(counts_set) = counts$X1
# Crea la tabla de información de las muestras targets_set con únicamente 
# las filas correspondientes a las muestras seleccionadas en counts_set.
targets_set = subset(targets, (targets$Sample_Name %in% df))
# Ordena alfanuméricamente las filas de targets_set en función de los 
# nombres de la columna Sample_Name para que éstas se correspondan con 
# el orden de las columnas de counts_set.
targets_set = targets_set[order(targets_set$Sample_Name),]


## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
library(DESeq2)
dds = DESeqDataSetFromMatrix(countData = counts_set, colData = targets_set,
                             design = ~ Group)


## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
nrow(dds)
dds = dds[rowSums(counts(dds)) > 1,]
nrow(dds)


## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
vsd = vst(dds, blind = FALSE)


## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
assay_vsd = t(assay(vsd))
rownames(assay_vsd) = vsd$Group
colnames(assay_vsd) = NULL
assay_vsd = assay_vsd[order(row.names(assay_vsd)),]
sampleDists = dist(assay_vsd)
library(pheatmap)
sampleDistMatrix = as.matrix(sampleDists)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         fontsize_row = 8, fontsize_col = 8,
         cluster_rows = T, cluster_cols = T)


## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
plotPCA(vsd, intgroup = "Group")


## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
library(ggplot2)
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(vsd)))
ggplot(mds, aes(X1,X2,color=Group)) + geom_point(size=3)


## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
dds = DESeq(dds, parallel =TRUE)

## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
res_SFI_NIT = results(dds, contrast = c("Group", "SFI", "NIT"), alpha = 0.05,
                      lfcThreshold = log2(2))
res_ELI_NIT = results(dds, contrast = c("Group", "ELI", "NIT"), alpha = 0.05,
                      lfcThreshold = log2(2))
res_ELI_SFI = results(dds, contrast = c("Group", "ELI", "SFI"), alpha = 0.05,
                      lfcThreshold = log2(2))


## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
summary(res_SFI_NIT)


## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
summary(res_ELI_NIT)


## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
summary(res_ELI_SFI)


## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
topGene_SFI_NIT = rownames(res_SFI_NIT)[which.min(res_SFI_NIT$padj)]
plotCounts(dds, topGene_SFI_NIT, "Group", col = "red3", pch = 16)

topGene_ELI_NIT = rownames(res_ELI_NIT)[which.min(res_ELI_NIT$padj)]
plotCounts(dds, topGene_ELI_NIT, "Group", col = "red3", pch = 16)

topGene_ELI_SFI = rownames(res_ELI_SFI)[which.min(res_ELI_SFI$padj)]
plotCounts(dds, topGene_ELI_SFI, "Group", col = "red3", pch = 16)

## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
plotMA(res_SFI_NIT, ylim=c(-25,25), colSig = "red3",
       colLine = rgb(1,0,0,.5))

plotMA(res_ELI_NIT, ylim=c(-25,25), colSig = "red3",
       colLine = rgb(1,0,0,.5))

plotMA(res_ELI_SFI, ylim=c(-10,10), colSig = "red3",
       colLine = rgb(1,0,0,.5))


## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
library("AnnotationDbi")
library(EnsDb.Hsapiens.v86)

res_SFI_NIT$symbol = mapIds(EnsDb.Hsapiens.v86,
                            keys=row.names(res_SFI_NIT),
                            column="SYMBOL",
                            keytype="GENEID",
                            multiVals="first")
res_SFI_NIT$entrezid = mapIds(EnsDb.Hsapiens.v86,
                            keys=row.names(res_SFI_NIT),
                            column="ENTREZID",
                            keytype="GENEID",
                            multiVals="first")

res_ELI_NIT$symbol = mapIds(EnsDb.Hsapiens.v86,
                            keys=row.names(res_ELI_NIT),
                            column="SYMBOL",
                            keytype="GENEID",
                            multiVals="first")
res_ELI_NIT$entrezid = mapIds(EnsDb.Hsapiens.v86,
                            keys=row.names(res_ELI_NIT),
                            column="ENTREZID",
                            keytype="GENEID",
                            multiVals="first")

res_ELI_SFI$symbol = mapIds(EnsDb.Hsapiens.v86,
                            keys=row.names(res_ELI_SFI),
                            column="SYMBOL",
                            keytype="GENEID",
                            multiVals="first")
res_ELI_SFI$entrezid = mapIds(EnsDb.Hsapiens.v86,
                            keys=row.names(res_ELI_SFI),
                            column="ENTREZID",
                            keytype="GENEID",
                            multiVals="first")


## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
topVarGenes = head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)

mat  = assay(vsd)[topVarGenes, ]
mat  = mat - rowMeans(mat)
colnames(mat) = targets_set$ShortName
row.names(mat) = mapIds(EnsDb.Hsapiens.v86,
                        keys=row.names(mat),
                        column="SYMBOL",
                        keytype="GENEID",
                        multiVals="first")
anno = as.data.frame(colData(vsd)[,"Group"])
row.names(anno) = colnames(mat)
colnames(anno) = "Group"
pheatmap(mat, annotation_col = anno, cluster_rows = T, cluster_cols = T,
         fontsize_row = 8, fontsize_col = 8)


## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
res_SFI_NIT_df = as.data.frame(res_SFI_NIT)
res_SFI_NIT_df = subset(res_SFI_NIT_df, padj < 0.05
                        & abs(log2FoldChange) > 1)
res_SFI_NIT_df = res_SFI_NIT_df[order(res_SFI_NIT_df$padj),]

res_ELI_NIT_df = as.data.frame(res_ELI_NIT)
res_ELI_NIT_df = subset(res_ELI_NIT_df, padj < 0.05
                        & abs(log2FoldChange) > 1)
res_ELI_NIT_df = res_ELI_NIT_df[order(res_ELI_NIT_df$padj),]

res_ELI_SFI_df = as.data.frame(res_ELI_SFI)
res_ELI_SFI_df = subset(res_ELI_SFI_df, padj < 0.05
                        & abs(log2FoldChange) > 1)
res_ELI_SFI_df = res_ELI_SFI_df[order(res_ELI_SFI_df$padj),]


## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
write.csv(res_SFI_NIT_df, file = "resultados/res_SFI_NIT_df.csv",
          row.names = T)

write.csv(res_ELI_NIT_df, file = "resultados/res_ELI_NIT_df.csv",
          row.names = T)

write.csv(res_ELI_SFI_df, file = "resultados/res_ELI_SFI_df.csv",
          row.names = T)


## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
sum(is.na(res_SFI_NIT_df$symbol))
sum(is.na(res_SFI_NIT_df$entrezid))

sum(is.na(res_ELI_NIT_df$symbol))
sum(is.na(res_ELI_NIT_df$entrezid))

sum(is.na(res_ELI_SFI_df$symbol))
sum(is.na(res_ELI_SFI_df$entrezid))


## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
res_SFI_NIT_genes = row.names(res_SFI_NIT_df)
res_ELI_NIT_genes = row.names(res_ELI_NIT_df)
res_ELI_SFI_genes = row.names(res_ELI_SFI_df)

comb = c(res_SFI_NIT_genes, res_ELI_NIT_genes, res_ELI_SFI_genes)

res_SFI_NIT_genes2 = comb %in% res_SFI_NIT_genes
res_ELI_NIT_genes2 = comb %in% res_ELI_NIT_genes
res_ELI_SFI_genes2 = comb %in% res_ELI_SFI_genes

counts_vennDiagram = cbind(res_SFI_NIT_genes2, res_ELI_NIT_genes2,
                           res_ELI_SFI_genes2)

library(limma)

results_vennDiagram = vennCounts(counts_vennDiagram)

for (i in 1:nrow(results_vennDiagram)) {
  d = sum(results_vennDiagram[i,-ncol(results_vennDiagram)])
  if (d == 0) next
  results_vennDiagram[i,"Counts"] = results_vennDiagram[i,"Counts"] / d
}

vennDiagram(results_vennDiagram, cex = 0.9, main = "",
            names = c("SFI vs NIT", "ELI vs NIT", "ELI vs SFI"),
            circle.col = c("red", "blue", "green3"))


## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
res_SFI_NIT_df = as.data.frame(res_SFI_NIT)
res_SFI_NIT_df = subset(res_SFI_NIT_df, padj < 0.15)

res_ELI_NIT_df = as.data.frame(res_ELI_NIT)
res_ELI_NIT_df = subset(res_ELI_NIT_df, padj < 0.15)

res_ELI_SFI_df = as.data.frame(res_ELI_SFI)
res_ELI_SFI_df = subset(res_ELI_SFI_df, padj < 0.15)

library(clusterProfiler)

ego_ELI_NIT = enrichGO(res_ELI_NIT_df$entrezid,
                       OrgDb = "org.Hs.eg.db",
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable = TRUE)

ego_ELI_SFI = enrichGO(res_ELI_SFI_df$entrezid,
                       OrgDb = "org.Hs.eg.db",
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable = TRUE)





## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'--------------------
ego_ELI_NIT = as.data.frame(ego_ELI_NIT)
write.csv(ego_ELI_NIT, file = "./resultados/ego_ELI_NIT.csv",
          row.names = F)

ego_ELI_SFI = as.data.frame(ego_ELI_SFI)
write.csv(ego_ELI_SFI, file = "./resultados/ego_ELI_SFI.csv",
          row.names = F)


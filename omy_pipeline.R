#library(rtracklayer)
library(tidyverse)
library(tximport)
library(matrixStats)
library(DESeq2)
library(limma)
library(edgeR)
#library(gt)
#library(DT)
#library(plotly)


library(gplots)
library(RColorBrewer)

#library(ensembldb)


targets <- read_tsv('studydesign_embrions.txt')

path <- file.path('transcripts_quant', 'embrions', targets$sample, 'quant.sf')
all(file.exists(path))

#tx <- read_tsv('../all_annotated_transcripts.csv')
#tx <- as_tibble(tx)

library(biomaRt)

listMarts()
#choose the 'mart' you want to work with
myMart <- useMart(biomart="ENSEMBL_MART_ENSEMBL")
#take a look at all available datasets within the selected mart
available.datasets <- listDatasets(myMart)
#now grab the ensembl annotations for dog
mykiss.anno <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "omykiss_gene_ensembl")
mykiss.attributes <- listAttributes(mykiss.anno)

Tx <- getBM(attributes=c(
                        'ensembl_transcript_id_version',
  #                      'external_gene_name',
                         'ensembl_gene_id'),
                mart = mykiss.anno)

#write.csv(Tx, 'gene_name_translator.csv')


Tx <- as_tibble(Tx)
#we need to rename the two columns we just retreived from biomart
Tx <- dplyr::rename(Tx, target_id = ensembl_transcript_id_version, 
                        gene_name = ensembl_gene_id)



Txi_gene <- tximport(path, 
                     type = "salmon", 
                     tx2gene = Tx, 
                     txOut = FALSE, #How does the result change if this =FALSE vs =TRUE?
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = FALSE
                     #ignoreAfterBar = TRUE
)

myTPM <- Txi_gene$abundance
myCounts <- Txi_gene$counts
colSums(myTPM)
colSums(myCounts)



#---------------
library(readr)
myCounts.df <- as_tibble(myCounts, rownames='geneID')


colnames(myCounts.df) <- c('geneID', sampleLabels)

write_tsv(myCounts.df, 'myCounts_do_clusta.txt')

#----------------



sampleLabels <- targets$sample

myTPM.stats <- transform(myTPM, 
                         SD=rowSds(myTPM), 
                         AVG=rowMeans(myTPM),
                         MED=rowMedians(myTPM))

head(myTPM.stats)

myDGEList <- DGEList(myCounts)

cpm <- cpm(myDGEList) 
log2.cpm <- cpm(myDGEList, log=TRUE)
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")


colnames(log2.cpm.df) <- c('geneID', sampleLabels)

#-------------------------

log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = sampleLabels, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

#-----------------------------------------------



keepers <- rowSums(cpm>1)>=4
myDGEList.filtered <- myDGEList[keepers,]

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
myDGEList.filtered.norm.df <- as_tibble(myDGEList.filtered.norm$counts, rownames = "geneID")
colnames(myDGEList.filtered.norm.df) <- c('geneID', targets$sample)

write.csv(myDGEList.filtered.norm.df, 'TMM_eggs_only.csv')
#write_tsv(myDGEList.filtered.norm.df, 'TMM_embrions_only.tsv')

# PIVOT ZNORMALIZOWANEJ EKSPRESJI DO ANOVY
cpm.filtered.norm.df.pivot <- pivot_longer(myDGEList.filtered.norm.df, # dataframe to be pivoted
                                                cols = sampleLabels, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)


write.csv(cpm.filtered.norm.df.pivot, 'TMM_embrions_only_PIVOT_ANOVA.csv')
# KONIEC PIVOT DO ZNOR. EKSP>


log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)


log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = sampleLabels, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

# FILTERED NORMALIZED
p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()













log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")


pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)

pc.var<-pca.res$sdev^2
pc.per<-round(pc.var/sum(pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
pc.per

pca.res.df <- as_tibble(pca.res$x)

ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color=targets$group) +
  geom_point(size=4) +
  #geom_label() +
  # stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  # coord_fixed() +
  theme_bw()



group <- factor(targets$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)


v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.filtered.norm, design)
contrast.matrix <- makeContrasts(infection = '17' - '18',
                                 levels=design)



# DESEQ2

targets$condition <- factor(targets$condition)


ddsTxi <- DESeqDataSetFromTximport(Txi_gene,
                                   colData = targets,
                                   design = ~ condition)
dds <- DESeq(ddsTxi)
res <- results(dds)
res <- results(dds, contrast=c("condition","KxZ","ZxZ"))
summary(res)
counts(dds)

#plotMA(res, ylim=c(-10,25))
write.csv(as.data.frame(res), file="deseq2_KxZ_vs_ZxZ.csv")

write.csv(as.data.frame(counts(dds)), file='deseq2_eggs_counts.csv')

counts_dds <- read.csv('deseq2_eggs_counts.csv')
colnames(counts_dds) <- c('geneID', sampleLabels)
write.csv(counts_dds, "deseq2_eggs_counts.csv")



# LRT----

ddsTxi <- DESeqDataSetFromTximport(Txi_gene,
                                   colData = targets,
                                   design = ~ condition)


dds <- DESeq(ddsTxi)
res <- results(dds)
counts(dds)
summary(res)



dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)
res_LRT <- results(dds_lrt)

res_LRT_tb <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigLRT_genes <- res_LRT_tb %>% dplyr::filter(padj < 0.05)

keep <- sigLRT_genes$gene

TMM_embrions <- read.csv('do_przeslania/TMM_embrions_matryca.csv')

TMM_embrions_padj_selected <- TMM_embrions[TMM_embrions$gene %in% keep,]

#zapis powyzszego, aby usunac kolumny do heatmapy
write.csv(TMM_embrions_padj_selected, file='TMM_LRT_selected_do_heatmapy.csv')

TMM_embrions_padj_selected_testy <- read.csv('TMM_LRT_selected_do_heatmapy.csv')

TMM_embrions_padj_selected_testy_matrix <- as.matrix(TMM_embrions_padj_selected_testy)

clustRows <- hclust(as.dist(1-cor(t(TMM_embrions_padj_selected_testy), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(TMM_embrions_padj_selected_testy, method="spearman")), method="complete")

module.assign <- cutree(clustRows, k=4)

#dodanie numerów klastrów do LRT
#TMM_LRT_selected <- read.csv('TMM_LRT_selected_eggs.csv') 
TMM_LRT_selected_w_clusters <- transform(TMM_embrions_padj_selected, 
                                        cluster.No=module.assign)
  #a następnie dodanie info o pval z testu LRT
sigLRT_genes <- dplyr::rename(sigLRT_genes,
                              'geneID'='gene')
TMM_LRT_selected_w_clusters_w_padj <- inner_join(TMM_LRT_selected_w_clusters, 
                                                 sigLRT_genes, by='geneID')
write.csv(TMM_LRT_selected_w_clusters_w_padj, 'TMM_LRT_selected_embrions.csv')
#end


module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 

myheatcolors3 <- brewer.pal(name="RdBu", n=11)

heatmap.2(TMM_embrions_padj_selected_testy_matrix, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=rev(myheatcolors3), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1) 
dev.off()


# heatmapa ANOVy z MARKIEM padj < 020 ----

anova_mark <- read.csv('d_p_val_embryos.csv')
anova_mark_padj_020 <- anova_mark[anova_mark$p_adjusted < 0.2,]
keep_mark_padj_020 <- anova_mark_padj_020$gene
TMM_embrions_padj_selected_anova_mark_padj_020 <- TMM_embrions[TMM_embrions$gene %in% keep_mark_padj_020,]
TMM_embrions_padj_selected_anova_mark_padj_020 <- subset(TMM_embrions_padj_selected_anova_mark_padj_020, select = -c(X, gene))
TMM_embrions_padj_selected_anova_mark_padj_020 <- as.matrix(TMM_embrions_padj_selected_anova_mark_padj_020)


clustRows <- hclust(as.dist(1-cor(t(TMM_embrions_padj_selected_anova_mark_padj_020), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(TMM_embrions_padj_selected_anova_mark_padj_020, method="spearman")), method="complete")

module.assign <- cutree(clustRows, k=4)

module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 

myheatcolors3 <- brewer.pal(name="RdBu", n=11)

heatmap.2(TMM_embrions_padj_selected_anova_mark_padj_020, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=rev(myheatcolors3), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1) 
dev.off()


# heatmapa ANOVy z MARKIEM padj < 025 ----

anova_mark <- read.csv('d_p_val_embryos.csv')
anova_mark_padj_025 <- anova_mark[anova_mark$p_adjusted < 0.25,]
keep_mark_padj_025 <- anova_mark_padj_025$gene
TMM_embrions_padj_selected_anova_mark_padj_025 <- TMM_embrions[TMM_embrions$gene %in% keep_mark_padj_025,]
TMM_embrions_padj_selected_anova_mark_padj_025 <- subset(TMM_embrions_padj_selected_anova_mark_padj_025, select = -c(X, gene))
TMM_embrions_padj_selected_anova_mark_padj_025 <- as.matrix(TMM_embrions_padj_selected_anova_mark_padj_025)


clustRows <- hclust(as.dist(1-cor(t(TMM_embrions_padj_selected_anova_mark_padj_025), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(TMM_embrions_padj_selected_anova_mark_padj_025, method="spearman")), method="complete")

module.assign <- cutree(clustRows, k=4)

module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 

myheatcolors3 <- brewer.pal(name="RdBu", n=11)

heatmap.2(TMM_embrions_padj_selected_anova_mark_padj_025, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=rev(myheatcolors3), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1) 
dev.off()






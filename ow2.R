library(dplyr)
library(openai)
library(Seurat)
library(cowplot)
library(ggplot2)
library(patchwork)
library(GPTCelltype)
library(org.Gg.eg.db)
library(EnhancedVolcano)
library(clusterProfiler)
Sys.setenv(OPENAI_API_KEY = '')
set.seed(42)

sample <- 'Ovary W2'
abbr <- 'OW2'
org <- 'Ovary'
setwd('/share/home/w/scRNA_seq/h/ovary')

t1 <- Read10X(data.dir = "../../01Cellranger_count_Results_2wks/01Count_Results_2wks/Ovary_R51_1_output/outs/filtered_feature_bc_matrix/")
t2 <- Read10X(data.dir = "../../01Cellranger_count_Results_2wks/01Count_Results_2wks/Ovary_R51_2_output/outs/filtered_feature_bc_matrix/")
t3 <- Read10X(data.dir = "../../01Cellranger_count_Results_2wks/01Count_Results_2wks/Ovary_R51_3_output/outs/filtered_feature_bc_matrix/")
c1 <- Read10X(data.dir = "../../01Cellranger_count_Results_2wks/01Count_Results_2wks/Ovary_PBS_1_output/outs/filtered_feature_bc_matrix/")
c2 <- Read10X(data.dir = "../../01Cellranger_count_Results_2wks/01Count_Results_2wks/Ovary_PBS_2_output/outs/filtered_feature_bc_matrix/")
c3 <- Read10X(data.dir = "../../01Cellranger_count_Results_2wks/01Count_Results_2wks/Ovary_PBS_3_output/outs/filtered_feature_bc_matrix/")
t1 <- CreateSeuratObject(counts = t1, project = paste(sample,"R51 R1"), min.cells = 3, min.features = 200)
t2 <- CreateSeuratObject(counts = t2, project = paste(sample,"R51 R2"), min.cells = 3, min.features = 200)
t3 <- CreateSeuratObject(counts = t3, project = paste(sample,"R51 R3"), min.cells = 3, min.features = 200)
c1 <- CreateSeuratObject(counts = c1, project = paste(sample,"PBS R1"), min.cells = 3, min.features = 200)
c2 <- CreateSeuratObject(counts = c2, project = paste(sample,"PBS R2"), min.cells = 3, min.features = 200)
c3 <- CreateSeuratObject(counts = c3, project = paste(sample,"PBS R3"), min.cells = 3, min.features = 200)

t1$group <- paste(sample,"R51")
t2$group <- paste(sample,"R51")
t3$group <- paste(sample,"R51")
c1$group <- paste(sample,"PBS")
c2$group <- paste(sample,"PBS")
c3$group <- paste(sample,"PBS")
t1$sample <- "R51 R1"
t2$sample <- "R51 R2"
t3$sample <- "R51 R3"
c1$sample <- "PBS R1"
c2$sample <- "PBS R2"
c3$sample <- "PBS R3"


#=====mt=====
t1[["percent.mt"]] <- PercentageFeatureSet(t1, pattern = "^(ND1|ND2|ND3|ND4|ND4L|ND5|ND6|COX1|COX2|COX3|ATP6|ATP8|CYTB)$")
t2[["percent.mt"]] <- PercentageFeatureSet(t2, pattern = "^(ND1|ND2|ND3|ND4|ND4L|ND5|ND6|COX1|COX2|COX3|ATP6|ATP8|CYTB)$")
t3[["percent.mt"]] <- PercentageFeatureSet(t3, pattern = "^(ND1|ND2|ND3|ND4|ND4L|ND5|ND6|COX1|COX2|COX3|ATP6|ATP8|CYTB)$")
c1[["percent.mt"]] <- PercentageFeatureSet(c1, pattern = "^(ND1|ND2|ND3|ND4|ND4L|ND5|ND6|COX1|COX2|COX3|ATP6|ATP8|CYTB)$")
c2[["percent.mt"]] <- PercentageFeatureSet(c2, pattern = "^(ND1|ND2|ND3|ND4|ND4L|ND5|ND6|COX1|COX2|COX3|ATP6|ATP8|CYTB)$")
c3[["percent.mt"]] <- PercentageFeatureSet(c3, pattern = "^(ND1|ND2|ND3|ND4|ND4L|ND5|ND6|COX1|COX2|COX3|ATP6|ATP8|CYTB)$")
t1[["percent.hb"]] <- PercentageFeatureSet(t1, pattern = "^(HBBA|HBAD|HBA1)$")
t2[["percent.hb"]] <- PercentageFeatureSet(t2, pattern = "^(HBBA|HBAD|HBA1)$")
t3[["percent.hb"]] <- PercentageFeatureSet(t3, pattern = "^(HBBA|HBAD|HBA1)$")
c1[["percent.hb"]] <- PercentageFeatureSet(c1, pattern = "^(HBBA|HBAD|HBA1)$")
c2[["percent.hb"]] <- PercentageFeatureSet(c2, pattern = "^(HBBA|HBAD|HBA1)$")
c3[["percent.hb"]] <- PercentageFeatureSet(c3, pattern = "^(HBBA|HBAD|HBA1)$")

# QC
t1 <- subset(t1, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 15 & percent.hb < 5)
t2 <- subset(t2, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 15 & percent.hb < 5)
t3 <- subset(t3, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 15 & percent.hb < 5)
c1 <- subset(c1, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 15 & percent.hb < 5)
c2 <- subset(c2, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 15 & percent.hb < 5)
c3 <- subset(c3, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 15 & percent.hb < 5)

t <- merge(t1, y = list(t2,t3))
c <- merge(c1, y = list(c2,c3))

t <- NormalizeData(t, verbose= FALSE)
c <- NormalizeData(c, verbose= FALSE)
t <- FindVariableFeatures(t, selection.method = "vst", nfeatures = 2000)
c <- FindVariableFeatures(c, selection.method = "vst", nfeatures = 2000)
t <- ScaleData(t, features = rownames(t))
c <- ScaleData(c, features = rownames(c))
t <- RunPCA(t, npcs = 100, verbose = FALSE)
c <- RunPCA(c, npcs = 100, verbose = FALSE)
t <- FindNeighbors(t, dims = 1:50, verbose = FALSE, reduction = "pca")
c <- FindNeighbors(c, dims = 1:50, verbose = FALSE, reduction = "pca")
t <- FindClusters(t,  resolution=0.5, verbose = FALSE)
c <- FindClusters(c,  resolution=0.5, verbose = FALSE)
t <- RunUMAP(t, dims = 1:50, verbose = FALSE, reduction = "pca")
c <- RunUMAP(c, dims = 1:50, verbose = FALSE, reduction = "pca")

#=====DimPlot without Cell ANNO=====
pdf(paste(abbr,'_R51_grid.pdf',sep=''), width = 30, height = 10)
p1 <- DimPlot(t, reduction = "umap", group.by = "sample") + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10))) + labs(title = paste(sample,'R51'))
p2 <- DimPlot(t, reduction = "umap", label = TRUE) + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10))) + labs(title = paste(sample,'R51'))
plot_grid(p1, p2, labels = c("A", "B"))
dev.off()
pdf(paste(abbr,'_R51_paired.pdf',sep=''), width = 40, height = 10)
DimPlot(t, reduction = "umap", split.by = "sample") + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
dev.off()

#=====DimPlot without Cell ANNO=====
pdf(paste(abbr,'_PBS_grid.pdf',sep=''), width = 30, height = 10)
p1 <- DimPlot(c, reduction = "umap", group.by = "sample") + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10))) + labs(title = paste(sample,'PBS'))
p2 <- DimPlot(c, reduction = "umap", label = TRUE) + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10))) + labs(title = paste(sample,'PBS'))
plot_grid(p1, p2, labels = c("A", "B"))
dev.off()
pdf(paste(abbr,'_PBS_paired.pdf',sep=''), width = 40, height = 10)
DimPlot(c, reduction = "umap", split.by = "sample") + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
dev.off()

# SCTransform
#install.packages('BiocManager')
#BiocManager::install('glmGamPoi')
paired.combined <- list(t1,t2,t3,c1,c2,c3)
paired.combined <- lapply(paired.combined, SCTransform, verbose = FALSE)
features <- SelectIntegrationFeatures(object.list = paired.combined, nfeatures = 2000)
paired.combined <- PrepSCTIntegration(object.list = paired.combined, anchor.features = features)
paired.anchors <- FindIntegrationAnchors(object.list = paired.combined, normalization.method = "SCT", anchor.features = features, dims = 1:50)
paired.combined <- IntegrateData(anchorset = paired.anchors, dims = 1:50, normalization.method = "SCT")
DefaultAssay(paired.combined) <- "integrated"
all.genes <- rownames(paired.combined)
paired.combined <- ScaleData(paired.combined, features = all.genes)
paired.combined <- RunPCA(paired.combined, npcs = 100, verbose = FALSE)
paired.combined <- FindNeighbors(paired.combined, dims = 1:50, verbose = FALSE, reduction = "pca")
paired.combined <- FindClusters(paired.combined, resolution=0.5, verbose = FALSE)
paired.combined <- RunUMAP(paired.combined, dims = 1:50, verbose = FALSE, reduction = "pca")
#paired.combined <- JoinLayers(paired.combined)
saveRDS(paired.combined, file = paste(abbr,"_object.rds",sep=''))

#=====DimPlot without Cell ANNO=====
pdf(paste(abbr,'_split.pdf',sep=''), width = 20, height = 10)
p1 <- DimPlot(paired.combined, reduction = "umap", group.by = "group") + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10))) + labs(title = sample)
p2 <- DimPlot(paired.combined, reduction = "umap", label = TRUE) + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10))) + labs(title = sample)
plot_grid(p1, p2, labels = c("A", "B"))
dev.off()

pdf(paste(abbr,'_paired.pdf',sep=''), width = 20, height = 10)
DimPlot(paired.combined, reduction = "umap", split.by = "group") + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
dev.off()

print('Find All Markers')
Idents(paired.combined) <- "seurat_clusters"
paired.markers <- FindAllMarkers(paired.combined, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

saveRDS(paired.markers, file = paste(abbr,"_markers.rds",sep=''))
write.csv(paired.markers, file = paste(abbr,'_genes_exp.csv',sep=''))

res <- gptcelltype(paired.markers, tissuename = paste('Chicken',org), model = 'gpt-4')
print(res)

#paired.combined <- readRDS(paste(abbr,"_object.rds",sep=''))

if (!dir.exists(paste0(abbr,"_DEGs"))) {
  dir.create(paste0(abbr,"_DEGs"))
}

if (!dir.exists(paste0(abbr,"_Volcano"))) {
  dir.create(paste0(abbr,"_Volcano"))
}

if (!dir.exists(paste0(abbr,"_KEGG"))) {
  dir.create(paste0(abbr,"_KEGG"))
}

# DEGs
all_clusters <- levels(Idents(paired.combined))
cluster_DEG_list <- list()
for (cluster_id in all_clusters) {
        # cluster cells
        cells_in_cluster <- WhichCells(paired.combined, idents = cluster_id)
        sub_obj <- subset(paired.combined, cells = cells_in_cluster)
        Idents(sub_obj) <- "group"
        deg <- tryCatch({
                FindMarkers(sub_obj, ident.1 = paste(sample,"R51"), ident.2 = paste(sample,"PBS"), logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox", only.pos = FALSE, min.cells.group = 10)
        }, error = function(e) {
                message(paste("Failed on cluster", cluster_id))
                return(NULL)
        })
        # CSV
        if (!is.null(deg)) {
                cluster_DEG_list[[cluster_id]] <- deg
                write.csv(deg, file = paste0(abbr,"_DEGs/", abbr, "_DEGs_c", cluster_id, "_R51vsPBS.csv"))
        }
}

# Volcano
for (cluster_id in names(cluster_DEG_list)) {
        df <- cluster_DEG_list[[cluster_id]]
        sig_genes <- df[df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25,]
        p <- EnhancedVolcano(df, lab = rownames(df), x = "avg_log2FC", y = "p_val_adj", pCutoff = 0.05, FCcutoff = 0.25, pointSize = 2.0, labSize = 2.0, title = paste("DEGs in Cluster", cluster_id), subtitle = paste0("Wilcoxon Test: R51 vs PBS, total DEGs: ", nrow(sig_genes)))
        print(p)
        ggsave(paste0(abbr,"_Volcano/", abbr, "_volcano_c", cluster_id, ".pdf"), p, width = 10, height = 8)
}


# KEGG
gga_kegg_local <- read.table("../kegg/gga_kegg_term2gene.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
term2gene <- gga_kegg_local[, c("pathway", "gene")]
term2name <- read.csv("../kegg/gga_kegg_term2name.csv", stringsAsFactors = FALSE)
kegg_result_list <- list()


for (cluster_id in names(cluster_DEG_list)) {
        deg <- cluster_DEG_list[[cluster_id]]
        # Filter out DEGs
        deg_filtered <- deg %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
        # Transfer symbol -> ENTREZ ID
        gene_df <- bitr(rownames(deg_filtered), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Gg.eg.db)
        gene_df <- gene_df[!duplicated(gene_df$ENTREZID) & !is.na(gene_df$ENTREZID), ]
        # KEGG
        cat("DEGs nrow: cluster", cluster_id, "\n")
        print(nrow(gene_df))
        if (nrow(gene_df) == 0) {
                cat("No valid ENTREZ IDs for cluster", cluster_id, "\n")
                next
        }

        kegg <- enricher(gene = gene_df$ENTREZID, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, TERM2GENE = term2gene, TERM2NAME = term2name)

        cat("KEGG nrow: cluster", cluster_id, "\n")
        if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
                kegg_result_list[[cluster_id]] <- kegg
                # Save as CSV
                write.csv(as.data.frame(kegg), paste0(abbr, "_KEGG/", abbr, "_KEGG_c", cluster_id, ".csv"), row.names = FALSE)
                # DotPlot
                pdf(paste0(abbr, "_KEGG/", abbr, "_KEGG_dot_c", cluster_id, ".pdf"))
                print(dotplot(kegg, showCategory = 20) + ggtitle(paste("Cluster", cluster_id)))
                dev.off()
        } else {
                cat("No enrichment result for cluster", cluster_id, "\n")
                # Write placeholder empty file
        }
}


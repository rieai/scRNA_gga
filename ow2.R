library(dplyr)
library(tidyr)
library(openai)
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(GPTCelltype)
library(org.Gg.eg.db)
library(clusterProfiler)
library(EnhancedVolcano)

Sys.setenv(OPENAI_API_KEY = '')
set.seed(42)

case <- 'Ovary W2'
abbr <- 'OW2'
org <- 'Ovary'
setwd('/share/home/wangbaikui/scRNA_seq/huanglei/ovary/ow2')

t1 <- Read10X(data.dir = "../../../01Cellranger_count_Results_2wks/01Count_Results_2wks/Ovary_R51_1_output/outs/filtered_feature_bc_matrix/")
t2 <- Read10X(data.dir = "../../../01Cellranger_count_Results_2wks/01Count_Results_2wks/Ovary_R51_2_output/outs/filtered_feature_bc_matrix/")
t3 <- Read10X(data.dir = "../../../01Cellranger_count_Results_2wks/01Count_Results_2wks/Ovary_R51_3_output/outs/filtered_feature_bc_matrix/")
c1 <- Read10X(data.dir = "../../../01Cellranger_count_Results_2wks/01Count_Results_2wks/Ovary_PBS_1_output/outs/filtered_feature_bc_matrix/")
c2 <- Read10X(data.dir = "../../../01Cellranger_count_Results_2wks/01Count_Results_2wks/Ovary_PBS_2_output/outs/filtered_feature_bc_matrix/")
c3 <- Read10X(data.dir = "../../../01Cellranger_count_Results_2wks/01Count_Results_2wks/Ovary_PBS_3_output/outs/filtered_feature_bc_matrix/")

t1 <- CreateSeuratObject(counts = t1, project = paste(abbr,"R51 R1"), min.cells = 3, min.features = 200)
t2 <- CreateSeuratObject(counts = t2, project = paste(abbr,"R51 R2"), min.cells = 3, min.features = 200)
t3 <- CreateSeuratObject(counts = t3, project = paste(abbr,"R51 R3"), min.cells = 3, min.features = 200)
c1 <- CreateSeuratObject(counts = c1, project = paste(abbr,"PBS R1"), min.cells = 3, min.features = 200)
c2 <- CreateSeuratObject(counts = c2, project = paste(abbr,"PBS R2"), min.cells = 3, min.features = 200)
c3 <- CreateSeuratObject(counts = c3, project = paste(abbr,"PBS R3"), min.cells = 3, min.features = 200)

t1$group <- paste(abbr, "R51")
t2$group <- paste(abbr, "R51")
t3$group <- paste(abbr, "R51")
c1$group <- paste(abbr, "PBS")
c2$group <- paste(abbr, "PBS")
c3$group <- paste(abbr, "PBS")
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
pdf(paste0(abbr,'_R51_grid.pdf'), width = 30, height = 10)
p1 <- DimPlot(t, reduction = "umap", group.by = "sample") + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10))) + labs(title = paste(case,'R51'))
p2 <- DimPlot(t, reduction = "umap", label = TRUE) + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10))) + labs(title = paste(case,'R51'))
plot_grid(p1, p2, labels = c("A", "B"))
dev.off()
pdf(paste0(abbr,'_R51_paired.pdf'), width = 40, height = 10)
DimPlot(t, reduction = "umap", split.by = "sample") + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
dev.off()

#=====DimPlot without Cell ANNO=====
pdf(paste0(abbr,'_PBS_grid.pdf'), width = 30, height = 10)
p1 <- DimPlot(c, reduction = "umap", group.by = "sample") + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10))) + labs(title = paste(case,'PBS'))
p2 <- DimPlot(c, reduction = "umap", label = TRUE) + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10))) + labs(title = paste(case,'PBS'))
plot_grid(p1, p2, labels = c("A", "B"))
dev.off()
pdf(paste0(abbr,'_PBS_paired.pdf'), width = 40, height = 10)
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
saveRDS(paired.combined, file = paste0('a_', abbr,"_object.rds"))

#=====DimPlot without Cell ANNO=====
pdf(paste0(abbr,'_split.pdf'), width = 20, height = 10)
p1 <- DimPlot(paired.combined, reduction = "umap", group.by = "group") + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10))) + labs(title = case)
p2 <- DimPlot(paired.combined, reduction = "umap", label = TRUE) + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10))) + labs(title = case)
plot_grid(p1, p2, labels = c("A", "B"))
dev.off()

pdf(paste0(abbr,'_paired.pdf'), width = 20, height = 10)
DimPlot(paired.combined, reduction = "umap", split.by = "group") + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
dev.off()

print('Find All Markers')
Idents(paired.combined) <- "seurat_clusters"
paired.markers <- FindAllMarkers(paired.combined, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
saveRDS(paired.markers, file = paste0('b_', abbr,"_markers.rds"))
write.csv(paired.markers, file = paste0(abbr,'_genes_exp.csv'))
res <- gptcelltype(paired.markers, tissuename = paste('Chicken',org), model = 'gpt-4')
print(res)

#paired.combined <- readRDS(paste0('a_', abbr, '_object.rds'))
#paired.markers <- readRDS(paste0('b_', abbr, '_markers.rds'))
paired.markers$cluster <- as.character(paired.markers$cluster)

# Top 20
print('Top 20 markers')
top20_genes_per_cluster <- paired.markers %>%
        group_by(cluster) %>%
        top_n(n = 20, wt = avg_log2FC) %>%
        arrange(cluster, desc(avg_log2FC)) %>%
        summarise(markers = paste(gene, collapse = ",")) %>%
        ungroup()

# Cell count
meta <- paired.combined@meta.data
df_counts <- as.data.frame(table(meta$seurat_clusters, meta$sample))
colnames(df_counts) <- c("cluster", "sample", "n_cells")
df_cell_counts <- pivot_wider(df_counts, names_from = sample, values_from = n_cells, values_fill = 0)
df_cell_counts$Cluster <- as.character(df_cell_counts$cluster)
tsv <- left_join(top20_genes_per_cluster, df_cell_counts, by = "cluster")

write.table(tsv, file = paste0(abbr, "_top20_markers.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

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
                print(paste0('FindMarkers ', cluster_id))
                FindMarkers(sub_obj, ident.1 = paste(abbr,"R51"), ident.2 = paste(abbr,"PBS"), logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox", only.pos = FALSE, min.cells.group = 10)
        }, error = function(e) {
                message(paste("Can not FindMarkers: cluster", cluster_id))
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
        up_genes <- sig_genes[sig_genes$avg_log2FC > 0.25, ]
        down_genes <- sig_genes[sig_genes$avg_log2FC < -0.25, ]
        top20_up <- up_genes[order(-up_genes$avg_log2FC), ][1:min(20, nrow(up_genes)), ]
        top20_down <- down_genes[order(down_genes$avg_log2FC), ][1:min(20, nrow(down_genes)), ]
        select_genes <- unique(c(rownames(top20_up), rownames(top20_down)))

        p <- EnhancedVolcano(df, lab = rownames(df), selectLab = select_genes, x = "avg_log2FC", y = "p_val_adj", drawConnectors = TRUE, pCutoff = 0.05, FCcutoff = 0.25, pointSize = 2.0, labSize = 3, title = paste("DEGs in Cluster", cluster_id), subtitle = paste0("Wilcoxon Test: R51 vs PBS, total DEGs: ", nrow(sig_genes)))
        ggsave(paste0(abbr,"_Volcano/", abbr, "_volcano_c", cluster_id, ".pdf"), p, width = 10, height = 8)
}

# KEGG
gga_kegg_local <- read.table("../../kegg/gga_kegg_term2gene.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
term2gene <- gga_kegg_local[, c("pathway", "gene")]
term2name <- read.csv("../../kegg/gga_kegg_term2name.csv", stringsAsFactors = FALSE)
kegg_result_list <- list()

all_DEGs <- c()
for (cluster_id in names(cluster_DEG_list)) {
        deg <- cluster_DEG_list[[cluster_id]]
        # Filter out DEGs
        deg_filtered <- deg %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
        all_DEGs <- c(all_DEGs, rownames(deg_filtered))
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


toll <- c('AKT1', 'AKT3', 'CASP18', 'CASP8', 'CCL4', 'CCL5', 'CD14', 'CD40', 'CD80', 'CD86', 'CD88', 'CHUK', 'CTSK', 'FADD', 'FOS', 'IFA3L', 'IFNA3', 'IFNAL4', 'IFNAL5', 'IFNAL6', 'IFNAR1', 'IFNAR2', 'IFNKL1', 'IFNW1', 'IKBKB', 'IKBKE', 'IL12A', 'IL12B', 'IL1B', 'IL6', 'IL8L1', 'IL8L2', 'IRAK4', 'IRF5', 'IRF7', 'JAK1', 'JUN', 'LOC100857744', 'LOC100857947', 'LOC100858177', 'LOC101750560', 'LOC395551', 'LOC768614', 'LY96', 'MAP2K1', 'MAP2K2', 'MAP2K3', 'MAP2K4', 'MAP2K6', 'MAP3K7', 'MAP3K8', 'MAPK1', 'MAPK10', 'MAPK11', 'MAPK12', 'MAPK13', 'MAPK14', 'MAPK3', 'MAPK8', 'MAPK9', 'MYD88', 'NFKB1', 'NFKBIA', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'RAC1', 'RELA', 'RIPK1', 'SPP1', 'STAT1', 'STAT2', 'TAB1', 'TAB2', 'TBK1', 'TICAM1', 'TIRAP', 'TLR1A', 'TLR1B', 'TLR2', 'TLR2A', 'TLR2B', 'TLR3', 'TLR4', 'TLR5', 'TLR7', 'TOLLIP', 'TRAF3', 'TRAF6', 'TYK2')
nod <- c('ANTXR1', 'ANTXR2', 'ANTXRL', 'ATG12', 'ATG16L1', 'ATG16L2', 'ATG5', 'AvBD2', 'BCL2', 'BCL2L1', 'BIRC2', 'BIRC8', 'BRCC3', 'CAMP', 'CARD8', 'CARD9', 'CASP1', 'CASP18', 'CASP8', 'CASR', 'CATH1', 'CATH2', 'CATH3', 'CCL5', 'CHUK', 'CTSB', 'CYBA', 'CYBB', 'DEFB4A', 'DHX33', 'DNM1L', 'ERBIN', 'FADD', 'GABARAPL1', 'GABARAPL2', 'GBP', 'GBP1', 'GBP4L', 'GBP7', 'GPRC6A', 'HSP90AA1', 'HSP90AB1', 'IFA3L', 'IFNA3', 'IFNAL4', 'IFNAL5', 'IFNAL6', 'IFNAR1', 'IFNAR2', 'IFNKL1', 'IFNW1', 'IKBKB', 'IKBKE', 'IL18', 'IL1B', 'IL6', 'IL8L1', 'IL8L2', 'IRAK4', 'IRF7', 'ITPR1', 'ITPR2', 'ITPR3', 'JAK1', 'JUN', 'LOC100857744', 'LOC100857947', 'LOC100858177', 'LOC107051192', 'LOC768614', 'MAP1LC3A', 'MAP1LC3B', 'MAP1LC3B2', 'MAP1LC3C', 'MAP3K7', 'MAPK1', 'MAPK10', 'MAPK11', 'MAPK12', 'MAPK13', 'MAPK14', 'MAPK3', 'MAPK8', 'MAPK9', 'MAVS', 'MCU', 'MFN1', 'MFN2', 'MYD88', 'NAMPT', 'NAMPTP1', 'NEK7', 'NFKB1', 'NFKBIA', 'NFKBIB', 'NLRP3', 'NLRPL', 'NLRX1', 'NOD1', 'P2RX7', 'PANX1', 'PKN2', 'PLCB1', 'PLCB2', 'PLCB4', 'PRKCD', 'PSTPIP1', 'RELA', 'RHOA', 'RIPK1', 'RIPK2', 'RIPK3', 'RNASEL', 'STAT1', 'STAT2', 'SUGT1', 'TAB1', 'TAB2', 'TAB3', 'TANK', 'TBK1', 'TICAM1', 'TLR4', 'TMEM173', 'TNFAIP3', 'TP53BP1', 'TRAF2', 'TRAF3', 'TRAF5', 'TRAF6', 'TRPM2', 'TRPM7', 'TRPV2', 'TXN', 'TXN2', 'TXNIP', 'TYK2', 'VDAC1', 'VDAC2', 'VDAC3', 'XIAP', 'YWHAE')
rig <- c('ADAR', 'ATG12', 'ATG5', 'AZI2', 'CASP10', 'CASP18', 'CASP8', 'CHUK', 'CYLD', 'DDX3X', 'DHX58', 'FADD', 'IFA3L', 'IFIH1', 'IFNA3', 'IFNAL4', 'IFNAL5', 'IFNAL6', 'IFNKL1', 'IFNW1', 'IKBKB', 'IKBKE', 'IL12A', 'IL12B', 'IL8L1', 'IL8L2', 'IRF7', 'LOC100857744', 'LOC100857947', 'LOC100858177', 'LOC100858279', 'LOC420294', 'LOC768614', 'MAP3K1', 'MAP3K7', 'MAPK10', 'MAPK11', 'MAPK12', 'MAPK13', 'MAPK14', 'MAPK8', 'MAPK9', 'MAVS', 'NFKB1', 'NFKBIA', 'NFKBIB', 'NLRX1', 'OTUD5', 'PIN1', 'RELA', 'RIPK1', 'SIKE1', 'TANK', 'TBK1', 'TBKBP1', 'TKFC', 'TMEM173', 'TRADD', 'TRAF2', 'TRAF3', 'TRAF6', 'TRIM25', 'ZNFX1')
cyto <- c('ADAR', 'CASP1', 'CASP18', 'CASP3', 'CASP7', 'CASP8', 'CCL4', 'CCL5', 'CGAS', 'CHUK', 'DDX41', 'DFNA5', 'DNASE2B', 'FADD', 'G3BP1', 'GSDME', 'IFA3L', 'IFNA3', 'IFNAL4', 'IFNAL5', 'IFNAL6', 'IFNKL1', 'IFNW1', 'IKBKB', 'IKBKE', 'IL18', 'IL1B', 'IL6', 'IRF7', 'LOC100857744', 'LOC100857947', 'LOC100858177', 'LOC395551', 'LOC768614', 'MAVS', 'MB21D1', 'MLKL', 'NFKB1', 'NFKBIA', 'NFKBIB', 'POLR1C', 'POLR1D', 'POLR2E', 'POLR2F', 'POLR2H', 'POLR2K', 'POLR2L', 'POLR3A', 'POLR3B', 'POLR3D', 'POLR3E', 'POLR3F', 'POLR3G', 'POLR3H', 'RELA', 'RIPK1', 'SAMHD1', 'TBK1', 'TMEM173', 'ZDHHC1')
lectin <- c('AKT1', 'AKT3', 'ARHGEF12', 'BCL10', 'CALM1', 'CALM2', 'CALML3', 'CALML4', 'CARD9', 'CASP1', 'CASP18', 'CASP8', 'CBLB', 'CCL17', 'CHUK', 'CYLD', 'EGR3', 'FCER1G', 'HRAS', 'IKBKB', 'IKBKE', 'IL10', 'IL12A', 'IL12B', 'IL1B', 'IL2', 'IL23A', 'IL6', 'IRF1', 'ITPR1', 'ITPR2', 'ITPR3', 'JUN', 'KRAS', 'KSR1', 'LOC101748577', 'LOC101748851', 'LOC101750560', 'LOC396477', 'LOC419429', 'LOC420294', 'LSP1', 'LSP1P1', 'MAP3K14', 'MAPK1', 'MAPK10', 'MAPK11', 'MAPK12', 'MAPK13', 'MAPK14', 'MAPK3', 'MAPK8', 'MAPK9', 'MAPKAPK2', 'MDM2', 'MRAS', 'NFATC1', 'NFATC2', 'NFATC3', 'NFKB1', 'NFKB2', 'NFKBIA', 'NRAS', 'PAK1', 'PCASP1', 'PCASP2', 'PCASP3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PLCG2', 'PLK3', 'PPP3CA', 'PPP3CB', 'PPP3CC', 'PPP3R1', 'PRKCD', 'PTGS2', 'PTPN11', 'RAF1', 'RELA', 'RHOA', 'RRAS', 'RRAS2', 'SRC', 'STAT1', 'STAT2', 'SYK')
iga <- c('AICDA', 'BLA', 'BLB1', 'BLB2', 'BLB3', 'CCR10', 'CCR9', 'CD28', 'CD40', 'CD40LG', 'CD80', 'CD86', 'CD88', 'CXCL12', 'CXCR4', 'DMA', 'DMB2', 'ICOS', 'ICOSLG', 'IL10', 'IL15', 'IL15RA', 'IL2', 'IL4', 'IL5', 'IL6', 'ITGA4', 'LOC107049666', 'MADCAM1', 'MAP3K14', 'MHCBL2', 'MHCDMA', 'MHCY2B1', 'MHCY2B2', 'MHCY2B6P', 'PIGR', 'TGFB1', 'TNFRSF13B', 'TNFRSF13C', 'TNFSF13B')
#05132  Salmonella infection
#05100  Bacterial invasion of epithelial cells
sal <- c('ABI1', 'ACBD3', 'ACTB', 'ACTG1', 'ACTG1L', 'ACTR10', 'ACTR10L', 'ACTR1A', 'ACTR2', 'ACTR3', 'ACTR3B', 'AHNAK', 'AHNAK2', 'AKT1', 'AKT3', 'ANXA2', 'ARF1', 'ARF6', 'ARHGEF26', 'ARL8A', 'ARL8B', 'ARL8BL', 'ARPC1A', 'ARPC1B', 'ARPC2', 'ARPC3', 'ARPC4', 'ARPC5', 'ARPC5L', 'BAK1', 'BCL2', 'BIRC2', 'BRK1', 'CASP1', 'CASP18', 'CASP3', 'CASP7', 'CASP8', 'CD14', 'CDC42', 'CHUK', 'CSE1L', 'CTNNB1', 'CYCS', 'CYFIP1', 'CYFIP2', 'CYTH1', 'CYTH3', 'CYTH4', 'DCTN1', 'DCTN2', 'DCTN3', 'DCTN4', 'DCTN5', 'DCTN6', 'DNM2L', 'DYL1', 'DYL2', 'DYNC1H1', 'DYNC1I1', 'DYNC1I2', 'DYNC1LI1', 'DYNC1LI2', 'DYNC1LI2L', 'DYNC2H1', 'DYNC2LI1', 'DYNLL1', 'DYNLL2', 'DYNLRB1', 'DYNLRB2', 'DYNLT1', 'DYNLT3', 'ELMO1', 'ELMO2', 'EXOC2', 'EXOC4', 'EXOC5', 'EXOC7', 'FADD', 'FBXO22', 'FHOD1', 'FLNA', 'FOS', 'FYCO1', 'GAPDH', 'GCC2', 'HRAS', 'HSP90AA1', 'HSP90AB1', 'HSP90B1', 'IKBKB', 'IL18', 'IL1B', 'IL6', 'IL8L1', 'IL8L2', 'IRAK4', 'JUN', 'KIF5B', 'KIF5C', 'KLC1', 'KLC2', 'KLC4', 'KPNA1', 'KPNA3', 'KPNA4', 'LEF1', 'LOC100857858', 'LOC100858984', 'LOC100859737', 'LOC101750560', 'LOC107050879', 'LOC107053365', 'LOC416695', 'LOC776720', 'LY96', 'M6PR', 'MAP2K1', 'MAP2K2', 'MAP2K3', 'MAP2K4', 'MAP2K6', 'MAP3K7', 'MAPK1', 'MAPK10', 'MAPK11', 'MAPK12', 'MAPK13', 'MAPK14', 'MAPK3', 'MAPK8', 'MAPK9', 'MLKL', 'MYC', 'MYD88', 'MYL10', 'MYL12A', 'MYL12B', 'MYL2', 'MYL9', 'MYLPF', 'MYO6', 'NCKAP1', 'NCKAP1L', 'NDAHNAKL', 'NFKB1', 'NFKBIA', 'NOD1', 'PAK1', 'PAK3', 'PFN1', 'PFN2', 'PFN3', 'PFN4', 'PIK3C2A', 'PIK3C2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PLEKHM1', 'PLEKHM2', 'PODXL', 'PTPRC', 'RAB5A', 'RAB5B', 'RAB5C', 'RAB7A', 'RAB7B', 'RAB9A', 'RAB9B', 'RAC1', 'RAF1', 'RALA', 'RELA', 'RHOA', 'RHOB', 'RHOG', 'RHOG2', 'RHOGL', 'RHOH', 'RHOJ', 'RILP', 'RIPK1', 'RIPK2', 'RIPK3', 'ROCK2', 'RPS3', 'S100A10', 'SKP1', 'SNX33', 'SNX9', 'TAB1', 'TAB2', 'TAB3', 'TCF7', 'TCF7L1', 'TCF7L2', 'TIRAP', 'TLR1A', 'TLR2', 'TLR2A', 'TLR2B', 'TLR4', 'TLR5', 'TNFRSF10B', 'TNFRSF1A', 'TNFSF10', 'TRADD', 'TRAF2', 'TRAF6', 'TUB4A', 'TUB5A', 'TUBA1A', 'TUBA1A1', 'TUBA1C', 'TUBA3E', 'TUBA4A', 'TUBA4AL', 'TUBA4B', 'TUBA8A', 'TUBA8B', 'TUBAL3', 'TUBB', 'TUBB1', 'TUBB2A', 'TUBB2B', 'TUBB3', 'TUBB4B', 'TUBB6', 'TXN', 'TXN2', 'VPS11', 'VPS16', 'VPS18', 'VPS33A', 'VPS39', 'VPS41', 'WASF3', 'WASF3L', 'WASL')
ec <- c('ACTB', 'ACTG1', 'ACTG1L', 'ACTR2', 'ACTR3', 'ACTR3B', 'ARHGAP10', 'ARHGEF26', 'ARPC1A', 'ARPC1B', 'ARPC2', 'ARPC3', 'ARPC4', 'ARPC5', 'ARPC5L', 'BCAR1', 'CAV1', 'CAV2', 'CAV3', 'CBL', 'CD2AP', 'CDC42', 'CDH1', 'CLTA', 'CLTB', 'CLTC', 'CLTCL1', 'CRK', 'CRKL', 'CTNNA1', 'CTNNA2', 'CTNNA3', 'CTNNB1', 'CTTN', 'DNM1', 'DNM2L', 'DNM3', 'DOCK1', 'ELMO1', 'ELMO2', 'ELMO3', 'FN1', 'GAB1', 'HCLS1', 'ILK', 'ITGB1', 'MAD2L2', 'MET', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PTK2', 'PXN', 'RAC1', 'RHOA', 'RHOG', 'RHOG2', 'RHOGL', 'SEPT11', 'SEPT12', 'SEPT2', 'SEPT2L', 'SEPT3', 'SEPT6', 'SEPT8', 'SEPT9', 'SEPTIN11', 'SEPTIN12', 'SEPTIN2', 'SEPTIN2L', 'SEPTIN3', 'SEPTIN6', 'SEPTIN8', 'SEPTIN9', 'SHC1', 'SHC2', 'SHC3', 'SHC4', 'SRC', 'VCL', 'WASF1', 'WASF2', 'WASL')

all_DEGs <- unique(all_DEGs)

if (!dir.exists(paste0(abbr,"_Dot"))) {
  dir.create(paste0(abbr,"_Dot"))
}


paired.combined$group_cluster <- paste(paired.combined$seurat_clusters, paired.combined$group, sep = "_")

# Dot
markers.to.plot <- intersect(all_DEGs, toll)
p <- DotPlot(paired.combined, features = rev(markers.to.plot), cols = c("#BB2A2D", "#1F4C98"), dot.scale = 8, split.by = "group") + RotatedAxis()
p <- p + xlab("Toll-like receptor signaling pathway") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"), axis.text.y = element_text(size = 10, face = "bold"), legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 9), panel.grid = element_blank(), legend.position = "right")
ggsave(paste0(abbr, "_Dot/a_", abbr, "_Toll.pdf"), plot = p, width = 12, height = 15, dpi = 1200)

cells.use <- paired.combined@meta.data %>%
        tibble::rownames_to_column("cell_id") %>%
        group_by(group_cluster) %>%
        sample_n(size = min(50, n()), replace = FALSE) %>%
        pull(cell_id)
cell_heat <- subset(paired.combined, cells = cells.use)

# Heatmap
p_heat <- DoHeatmap(cell_heat, features = markers.to.plot, group.by = "group_cluster", size = 7) + scale_fill_gradientn(colors = c("#1F4C98", "white", "#BB2A2D")) + theme(axis.text.x = element_text(size = 5)) + ggtitle("Toll-like receptor signaling pathway")

ggsave(paste0(abbr,"_Dot/b_", abbr, "_Toll_Heat.pdf"), plot = p_heat, width = 20, height = 12, dpi = 1200)

# Violin
pdf(paste0(abbr,"_Dot/c_", abbr, "_Toll_vln.pdf"), width = 15, height = 5*length(markers.to.plot))
plots <- VlnPlot(paired.combined, features = markers.to.plot, split.by = "group", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()

# Feature
pdf(paste0(abbr,"_Dot/d_", abbr, "_Toll_feature.pdf"), width = 20, height = 5*length(markers.to.plot))
FeaturePlot(paired.combined, features = markers.to.plot, split.by = "group", cols = c("grey","red"), reduction = "umap", min.cutoff = "q9") + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
dev.off()

# Dot
markers.to.plot <- intersect(all_DEGs, nod)
p <- DotPlot(paired.combined, features = rev(markers.to.plot), cols = c("#BB2A2D", "#1F4C98"), dot.scale = 8, split.by = "group") + RotatedAxis()
p <- p + xlab("NOD-like receptor signaling pathway") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"), axis.text.y = element_text(size = 10, face = "bold"), legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 9), panel.grid = element_blank(), legend.position = "right")
ggsave(paste0(abbr, "_Dot/a_", abbr, "_NOD.pdf"), plot = p, width = 12, height = 15, dpi = 1200)

# Heatmap
p_heat <- DoHeatmap(cell_heat, features = markers.to.plot, group.by = "group_cluster", size = 6) + scale_fill_gradientn(colors = c("#1F4C98", "white", "#BB2A2D")) + theme(axis.text.x = element_text(size = 4)) + ggtitle("NOD-like receptor signaling pathway")
ggsave(paste0(abbr,"_Dot/b_", abbr, "_NOD_Heat.pdf"), plot = p_heat, width = 20, height = 12, dpi = 1200)


# Violin
pdf(paste0(abbr,"_Dot/c_", abbr, "_NOD_vln.pdf"), width = 15, height = 5*length(markers.to.plot))
plots <- VlnPlot(paired.combined, features = markers.to.plot, split.by = "group", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()

# Feature
pdf(paste0(abbr,"_Dot/d_", abbr, "_NOD_feature.pdf"), width = 20, height = 5*length(markers.to.plot))
FeaturePlot(paired.combined, features = markers.to.plot, split.by = "group", cols = c("grey","red"), reduction = "umap", min.cutoff = "q9") + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
dev.off()

# RIG-I-like receptor signaling pathway
markers.to.plot <- intersect(all_DEGs, rig)
p <- DotPlot(paired.combined, features = rev(markers.to.plot), cols = c("#BB2A2D", "#1F4C98"), dot.scale = 8, split.by = "group") + RotatedAxis()
p <- p + xlab("RIG-I-like receptor signaling pathway") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"), axis.text.y = element_text(size = 10, face = "bold"), legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 9), panel.grid = element_blank(), legend.position = "right")
ggsave(paste0(abbr, "_Dot/a_", abbr, "_RIG.pdf"), plot = p, width = 12, height = 15, dpi = 1200)

# Heatmap
p_heat <- DoHeatmap(cell_heat, features = markers.to.plot, group.by = "group_cluster", size = 6) + scale_fill_gradientn(colors = c("#1F4C98", "white", "#BB2A2D")) + theme(axis.text.x = element_text(size = 4)) + ggtitle("RIG-I-like receptor signaling pathway")

ggsave(paste0(abbr,"_Dot/b_", abbr, "_RIG_Heat.pdf"), plot = p_heat, width = 20, height = 12, dpi = 1200)


# Violin
pdf(paste0(abbr,"_Dot/c_", abbr, "_RIG_vln.pdf"), width = 15, height = 5*length(markers.to.plot))
plots <- VlnPlot(paired.combined, features = markers.to.plot, split.by = "group", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()

# Feature
pdf(paste0(abbr,"_Dot/d_", abbr, "_RIG_feature.pdf"), width = 20, height = 5*length(markers.to.plot))
FeaturePlot(paired.combined, features = markers.to.plot, split.by = "group", cols = c("grey","red"), reduction = "umap", min.cutoff = "q9") + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
dev.off()

# Cytosolic DNA-sensing pathway
markers.to.plot <- intersect(all_DEGs, cyto)
p <- DotPlot(paired.combined, features = rev(markers.to.plot), cols = c("#BB2A2D", "#1F4C98"), dot.scale = 8, split.by = "group") + RotatedAxis()
p <- p + xlab("Cytosolic DNA-sensing pathway") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"), axis.text.y = element_text(size = 10, face = "bold"), legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 9), panel.grid = element_blank(), legend.position = "right")
ggsave(paste0(abbr,"_Dot/a_", abbr, "_Cytosolic.pdf"), plot = p, width = 12, height = 15, dpi = 1200)

# Heatmap
p_heat <- DoHeatmap(cell_heat, features = markers.to.plot, group.by = "group_cluster", size = 6) + scale_fill_gradientn(colors = c("#1F4C98", "white", "#BB2A2D")) + theme(axis.text.x = element_text(size = 4)) + ggtitle("Cytosolic DNA-sensing pathway")

ggsave(paste0(abbr,"_Dot/b_", abbr, "_Cytosolic_Heat.pdf"), plot = p_heat, width = 20, height = 12, dpi = 1200)

# Violin
pdf(paste0(abbr,"_Dot/c_", abbr, "_Cytosolic_vln.pdf"), width = 15, height = 5*length(markers.to.plot))
plots <- VlnPlot(paired.combined, features = markers.to.plot, split.by = "group", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()

# Feature
pdf(paste0(abbr,"_Dot/d_", abbr, "_Cytosolic_feature.pdf"), width = 20, height = 5*length(markers.to.plot))
FeaturePlot(paired.combined, features = markers.to.plot, split.by = "group", cols = c("grey","red"), reduction = "umap", min.cutoff = "q9") + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
dev.off()

# C-type lectin receptor signaling pathway
markers.to.plot <- intersect(all_DEGs, lectin)
p <- DotPlot(paired.combined, features = rev(markers.to.plot), cols = c("#BB2A2D", "#1F4C98"), dot.scale = 8, split.by = "group") + RotatedAxis()
p <- p + xlab("C-type lectin receptor signaling pathway") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"), axis.text.y = element_text(size = 10, face = "bold"), legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 9), panel.grid = element_blank(), legend.position = "right")
ggsave(paste0(abbr,"_Dot/a_", abbr, "_Lectin.pdf"), plot = p, width = 12, height = 15, dpi = 1200)

# Heatmap
p_heat <- DoHeatmap(cell_heat, features = markers.to.plot, group.by = "group_cluster", size = 6) + scale_fill_gradientn(colors = c("#1F4C98", "white", "#BB2A2D")) + theme(axis.text.x = element_text(size = 4)) + ggtitle("C-type lectin receptor signaling pathway")
ggsave(paste0(abbr,"_Dot/b_", abbr, "_Lectin_Heat.pdf"), plot = p_heat, width = 20, height = 12, dpi = 1200)

# Violin
pdf(paste0(abbr,"_Dot/c_", abbr, "_Lectin_vln.pdf"), width = 15, height = 5*length(markers.to.plot))
plots <- VlnPlot(paired.combined, features = markers.to.plot, split.by = "group", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()

# Feature
pdf(paste0(abbr,"_Dot/d_", abbr, "_Lectin_feature.pdf"), width = 20, height = 5*length(markers.to.plot))
FeaturePlot(paired.combined, features = markers.to.plot, split.by = "group", cols = c("grey","red"), reduction = "umap", min.cutoff = "q9") + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
dev.off()

# Intestinal immune network for IgA production
markers.to.plot <- intersect(all_DEGs, iga)
p <- DotPlot(paired.combined, features = rev(markers.to.plot), cols = c("#BB2A2D", "#1F4C98"), dot.scale = 8, split.by = "group") + RotatedAxis()
p <- p + xlab("Intestinal immune network for IgA production") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"), axis.text.y = element_text(size = 10, face = "bold"), legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 9), panel.grid = element_blank(), legend.position = "right")
ggsave(paste0(abbr,"_Dot/a_", abbr, "_IGA.pdf"), plot = p, width = 12, height = 15, dpi = 1200)

# Heatmap
p_heat <- DoHeatmap(cell_heat, features = markers.to.plot, group.by = "group_cluster", size = 6) + scale_fill_gradientn(colors = c("#1F4C98", "white", "#BB2A2D")) + theme(axis.text.x = element_text(size = 4)) + ggtitle("Intestinal immune network for IgA production")
ggsave(paste0(abbr,"_Dot/b_", abbr, "_IGA_Heat.pdf"), plot = p_heat, width = 20, height = 12, dpi = 1200)

# Violin
pdf(paste0(abbr,"_Dot/c_", abbr, "_IGA_vln.pdf"), width = 15, height = 5*length(markers.to.plot))
plots <- VlnPlot(paired.combined, features = markers.to.plot, split.by = "group", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()

# Feature
pdf(paste0(abbr,"_Dot/d_", abbr, "_IGA_feature.pdf"), width = 20, height = 5*length(markers.to.plot))
FeaturePlot(paired.combined, features = markers.to.plot, split.by = "group", cols = c("grey","red"), reduction = "umap", min.cutoff = "q9") + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
dev.off()

# Salmonella infection
markers.to.plot <- intersect(all_DEGs, sal)
p <- DotPlot(paired.combined, features = rev(markers.to.plot), cols = c("#BB2A2D", "#1F4C98"), dot.scale = 8, split.by = "group") + RotatedAxis()
p <- p + xlab("Salmonella infection") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"), axis.text.y = element_text(size = 10, face = "bold"), legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 9), panel.grid = element_blank(), legend.position = "right")
ggsave(paste0(abbr,"_Dot/a_", abbr, "_Sal.pdf"), plot = p, width = 12, height = 15, dpi = 1200)

# Heatmap
p_heat <- DoHeatmap(cell_heat, features = markers.to.plot, group.by = "group_cluster", size = 6) + scale_fill_gradientn(colors = c("#1F4C98", "white", "#BB2A2D")) + theme(axis.text.x = element_text(size = 4)) + ggtitle("Salmonella infection")
ggsave(paste0(abbr,"_Dot/b_", abbr, "_Sal_Heat.pdf"), plot = p_heat, width = 20, height = 12, dpi = 1200)

# Violin
pdf(paste0(abbr,"_Dot/c_", abbr, "_Sal_vln.pdf"), width = 15, height = 5*length(markers.to.plot))
plots <- VlnPlot(paired.combined, features = markers.to.plot, split.by = "group", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()

# Feature
pdf(paste0(abbr,"_Dot/d_", abbr, "_Sal_feature.pdf"), width = 20, height = 5*length(markers.to.plot))
FeaturePlot(paired.combined, features = markers.to.plot, split.by = "group", cols = c("grey","red"), reduction = "umap", min.cutoff = "q9") + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
dev.off()

# Bacterial invasion of epithelial cells
markers.to.plot <- intersect(all_DEGs, ec)
p <- DotPlot(paired.combined, features = rev(markers.to.plot), cols = c("#BB2A2D", "#1F4C98"), dot.scale = 8, split.by = "group") + RotatedAxis()
p <- p + xlab("Bacterial invasion of epithelial cells") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"), axis.text.y = element_text(size = 10, face = "bold"), legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 9), panel.grid = element_blank(), legend.position = "right")
ggsave(paste0(abbr,"_Dot/a_", abbr, "_ECs.pdf"), plot = p, width = 12, height = 15, dpi = 1200)

# Heatmap
p_heat <- DoHeatmap(cell_heat, features = markers.to.plot, group.by = "group_cluster", size = 6) + scale_fill_gradientn(colors = c("#1F4C98", "white", "#BB2A2D")) + theme(axis.text.x = element_text(size = 4)) + ggtitle("Bacterial invasion of epithelial cells")
ggsave(paste0(abbr,"_Dot/b_", abbr, "_ECs_Heat.pdf"), plot = p_heat, width = 20, height = 12, dpi = 1200)

# Violin
pdf(paste0(abbr,"_Dot/c_", abbr, "_ECs_vln.pdf"), width = 15, height = 5*length(markers.to.plot))
plots <- VlnPlot(paired.combined, features = markers.to.plot, split.by = "group", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()

# Feature
pdf(paste0(abbr,"_Dot/d_", abbr, "_ECs_feature.pdf"), width = 20, height = 5*length(markers.to.plot))
FeaturePlot(paired.combined, features = markers.to.plot, split.by = "group", cols = c("grey","red"), reduction = "umap", min.cutoff = "q9") + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))

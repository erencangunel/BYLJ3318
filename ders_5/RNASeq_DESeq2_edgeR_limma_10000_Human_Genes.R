# ==============================================================================
# R İLE BİYOİNFORMATİK DERSİ UYGULAMA KODLARI (10.000 Gerçek İnsan Geni Seti)
# Konu: Ham RNA-Seq Count Verisinden DESeq2, edgeR ve limma-voom Analizlerine
# Model: İnsan Kolorektal Adenokarsinom Evreleri (Evre I, Evre II) ve FOLFOX Kemoterapi Yanıtı
# Örneklem: Homo_sapiens_genes.tsv İçinden Seçilen 10.000 Gerçek Protein Kodlayan İnsan Geni x 24 Örnek
# Özellikler:
#   * Her karşılaştırmada 180 ile 300 arası Up-regüle, 60 ile 120 arası Down-regüle Gen
#   * Karşılaştırmalar arası Up-regüle ve Down-regüle gen kümelerinde %10 ile %17 arası ortak (overlapped) genler
# ==============================================================================

library(DESeq2)
library(edgeR)
library(limma)
library(ggplot2)
library(pheatmap)

# ==============================================================================
# ADIM 1: Homo_sapiens_genes.tsv TABLOSUNDAN 10.000 GERÇEK İNSAN GENİNİN SEÇİMİ
#         VE COUNT MATRİSİNİN ÜRETİMİ (VEYA DOĞRUDAN CSV DOSYASINDAN OKUMA)
# ==============================================================================
set.seed(2026)

n_genes <- 10000
n_samples <- 24

group_levels <- c("Control", "Stage1", "Stage2", "PostChemo")
groups <- factor(rep(group_levels, each = 6), levels = group_levels)
sample_names <- paste0(rep(c("Ctrl", "Stg1", "Stg2", "Chemo"), each = 6), "_Rep", 1:6)

col_data <- data.frame(
  SampleID = sample_names,
  Condition = groups,
  PatientSex = factor(rep(c("M", "F"), 12)),
  row.names = sample_names
)

# Anahtar biyobelirteçlerin tanımlanması
known_oncogenes <- c("MYC", "KRAS", "CCND1", "CDK4", "PCNA", "MKI67")
known_tsg       <- c("APC", "TP53", "PTEN", "SMAD4", "CDH1", "RB1")
known_invasion  <- c("VEGFA", "MMP9", "HIF1A", "TWIST1", "SNAI1", "VIM")
known_chemo     <- c("CDKN1A", "BAX", "GADD45A", "TYMS", "FAS", "BBC3")
all_known       <- c(known_oncogenes, known_tsg, known_invasion, known_chemo)

# Eğer çalışma dizininde Homo_sapiens_genes.tsv varsa 10.000 gerçek insan gen sembolünü oradan çeker
if (file.exists("Homo_sapiens_genes.tsv")) {
  genes_tsv <- read.delim("Homo_sapiens_genes.tsv", stringsAsFactors = FALSE)
  valid_chroms <- c(as.character(1:22), "X", "Y")
  pc_genes <- subset(genes_tsv, type == "protein_coding" & chromosome %in% valid_chroms &
                     !is.na(gene_symbol) & gene_symbol != "")
  pc_genes <- pc_genes[!duplicated(pc_genes$gene_symbol), ]
  
  other_pool <- setdiff(pc_genes$gene_symbol, all_known)
  sampled_others <- sample(other_pool, size = n_genes - length(all_known), replace = FALSE)
  
  gene_names <- character(n_genes)
  gene_names[1:6]     <- known_oncogenes
  gene_names[226:231] <- known_invasion
  gene_names[461:466] <- known_chemo
  gene_names[646:651] <- known_tsg
  
  empty_idx <- which(gene_names == "")
  gene_names[empty_idx] <- sampled_others
} else if (file.exists("CRC_RNASeq_10000_Human_Genes_Counts.csv")) {
  existing_counts <- read.csv("CRC_RNASeq_10000_Human_Genes_Counts.csv", row.names = 1, check.names = FALSE)
  gene_names <- rownames(existing_counts)
} else {
  stop("Lütfen Homo_sapiens_genes.tsv veya CRC_RNASeq_10000_Human_Genes_Counts.csv dosyasını çalışma dizinine ekleyin.")
}

cat("Toplam seçilen gerçek insan geni sayısı:", length(gene_names), "\n")
cat("İlk 15 gen sembolü:", paste(head(gene_names, 15), collapse = ", "), "\n")

# Bazal ifade düzeyleri (İlk 8500 gen aktif, son 1500 gen düşük ifadeli gürültü)
base_means <- numeric(n_genes)
base_means[1:8500] <- round(rlnorm(8500, meanlog = 5.8, sdlog = 0.80))
base_means[1:8500][base_means[1:8500] < 45] <- 50
base_means[8501:10000] <- round(runif(1500, min = 0.5, max = 4.5))

fold_changes <- matrix(1, nrow = n_genes, ncol = 4)
colnames(fold_changes) <- group_levels

# Blok 1 (Gen 1 ile 36; 36 gen): Stage1_UP ve Stage2_UP Ortak Onkogenleri (~%13-%16 Overlap)
fold_changes[1:36, "Stage1"]    <- runif(36, 3.4, 4.8)
fold_changes[1:36, "Stage2"]    <- fold_changes[1:36, "Stage1"] * runif(36, 0.96, 1.05)
fold_changes[1:24, "PostChemo"] <- fold_changes[1:24, "Stage2"] * runif(24, 0.18, 0.26)
fold_changes[25:36, "PostChemo"] <- fold_changes[25:36, "Stage2"] * runif(12, 0.96, 1.04)

# Blok 2 (Gen 37 ile 225; 189 gen): Stage1 Spesifik UP (Toplam Stage1_UP = 225 gen; 28 tanesi PostChemo_UP ile ortak)
fold_changes[37:225, "Stage1"]    <- runif(189, 3.3, 4.7)
fold_changes[37:225, "Stage2"]    <- runif(189, 0.96, 1.04)
fold_changes[37:64, "PostChemo"]  <- runif(28, 3.5, 5.0)
fold_changes[65:225, "PostChemo"] <- runif(161, 0.96, 1.04)

# Blok 3 (Gen 226 ile 460; 235 gen): Stage2 Spesifik UP (Toplam Stage2_UP = 271 gen; 28 tanesi PostChemo_UP ile ortak)
fold_changes[226:460, "Stage1"]    <- runif(235, 0.96, 1.04)
fold_changes[226:460, "Stage2"]    <- runif(235, 3.3, 4.7)
fold_changes[226:273, "PostChemo"] <- fold_changes[226:273, "Stage2"] * runif(48, 0.18, 0.26)
fold_changes[274:301, "PostChemo"] <- fold_changes[274:301, "Stage2"] * runif(28, 3.3, 4.5)
fold_changes[302:460, "PostChemo"] <- fold_changes[302:460, "Stage2"] * runif(159, 0.96, 1.04)

# Blok 4 (Gen 461 ile 645; 185 gen): Saf PostChemo Spesifik UP (Toplam PostChemo_UP = 241 gen)
fold_changes[461:645, "Stage1"]    <- runif(185, 0.96, 1.04)
fold_changes[461:645, "Stage2"]    <- runif(185, 0.96, 1.04)
fold_changes[461:645, "PostChemo"] <- runif(185, 3.5, 5.0)

# Blok 5 (Gen 646 ile 659; 14 gen): Stage1_DOWN ve Stage2_DOWN Ortak Tümör Baskılayıcılar (~%13-%16 Overlap)
fold_changes[646:659, "Stage1"]    <- runif(14, 0.18, 0.27)
fold_changes[646:659, "Stage2"]    <- fold_changes[646:659, "Stage1"] * runif(14, 0.96, 1.04)
fold_changes[646:659, "PostChemo"] <- fold_changes[646:659, "Stage2"] * runif(14, 0.96, 1.04)

# Blok 6 (Gen 660 ile 730; 71 gen): Stage1 Spesifik DOWN (Toplam Stage1_DOWN = 85 gen; 11 tanesi PostChemo_DOWN ile ortak)
fold_changes[660:730, "Stage1"]    <- runif(71, 0.18, 0.27)
fold_changes[660:730, "Stage2"]    <- runif(71, 0.96, 1.04)
fold_changes[660:670, "PostChemo"] <- runif(11, 0.18, 0.26)
fold_changes[671:730, "PostChemo"] <- runif(60, 0.96, 1.04)

# Blok 7 (Gen 731 ile 826; 96 gen): Stage2 Spesifik DOWN (Toplam Stage2_DOWN = 110 gen; 12 tanesi PostChemo_DOWN ile ortak)
fold_changes[731:826, "Stage1"]    <- runif(96, 0.96, 1.04)
fold_changes[731:826, "Stage2"]    <- runif(96, 0.22, 0.30)
fold_changes[731:742, "PostChemo"] <- fold_changes[731:742, "Stage2"] * runif(12, 0.20, 0.26)
fold_changes[743:826, "PostChemo"] <- fold_changes[743:826, "Stage2"] * runif(84, 0.96, 1.04)

library_size_factors <- runif(n_samples, min = 0.85, max = 1.20)
sample_groups <- as.character(col_data$Condition)
mu_matrix <- (base_means * fold_changes[, sample_groups]) * rep(library_size_factors, each = n_genes)

counts_matrix <- matrix(
  rnbinom(n_genes * n_samples, mu = as.vector(mu_matrix), size = 1 / 0.07),
  nrow = n_genes,
  ncol = n_samples,
  dimnames = list(gene_names, sample_names)
)

# Matris boyutunun kontrolü (10.000 satır x 24 sütun)
cat("Count matrisi boyutu (Gen x Örnek):", dim(counts_matrix)[1], "x", dim(counts_matrix)[2], "\n")

write.csv(counts_matrix, "CRC_RNASeq_10000_Human_Genes_Counts.csv")
write.csv(col_data, "CRC_RNASeq_Metadata.csv")

# ==============================================================================
# ADIM 2: DESeq2 İLE ANALİZ VE GRUPLAR ARASI OVERLAP (%) HESAPLAMASI
# ==============================================================================
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = col_data,
  design = ~ Condition
)

keep <- rowSums(counts(dds) >= 10) >= 6
dds <- dds[keep, ]

cat("Filtreleme öncesi toplam gen sayısı:", nrow(counts_matrix), "\n")
cat("DESeq2 filtrelemesi sonrası kalan aktif gen sayısı:", nrow(dds), "\n")

dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "Condition") +
  ggtitle("PCA: Kolorektal Kanser Evreleri ve Kemoterapi Yanıtı (10.000 İnsan Geni)") +
  theme_bw()

res_Stg1_vs_Ctrl  <- as.data.frame(results(dds, contrast = c("Condition", "Stage1", "Control"), alpha = 0.05))
res_Stg2_vs_Ctrl  <- as.data.frame(results(dds, contrast = c("Condition", "Stage2", "Control"), alpha = 0.05))
res_Chemo_vs_Stg2 <- as.data.frame(results(dds, contrast = c("Condition", "PostChemo", "Stage2"), alpha = 0.05))

s1_up   <- rownames(subset(res_Stg1_vs_Ctrl,  padj < 0.05 & log2FoldChange > 1))
s1_down <- rownames(subset(res_Stg1_vs_Ctrl,  padj < 0.05 & log2FoldChange < -1))

s2_up   <- rownames(subset(res_Stg2_vs_Ctrl,  padj < 0.05 & log2FoldChange > 1))
s2_down <- rownames(subset(res_Stg2_vs_Ctrl,  padj < 0.05 & log2FoldChange < -1))

ch_up   <- rownames(subset(res_Chemo_vs_Stg2, padj < 0.05 & log2FoldChange > 1))
ch_down <- rownames(subset(res_Chemo_vs_Stg2, padj < 0.05 & log2FoldChange < -1))

cat(sprintf("Stage1 vs Control   -> UP: %d | DOWN: %d\n", length(s1_up), length(s1_down)))
cat(sprintf("Stage2 vs Control   -> UP: %d | DOWN: %d\n", length(s2_up), length(s2_down)))
cat(sprintf("PostChemo vs Stage2 -> UP: %d | DOWN: %d\n\n", length(ch_up), length(ch_down)))

check_overlap <- function(name1, set1, name2, set2) {
  common <- length(intersect(set1, set2))
  pct1 <- 100 * common / length(set1)
  pct2 <- 100 * common / length(set2)
  cat(sprintf("%s & %s -> %d ortak gen (%%%.1f / %%%.1f)\n", name1, name2, common, pct1, pct2))
}

cat("--- UP-REGÜLE GEN KÜMELERİ ARASI KESİŞİM (OVERLAP) ---\n")
check_overlap("Stage1_UP", s1_up, "Stage2_UP", s2_up)
check_overlap("Stage2_UP", s2_up, "PostChemo_UP", ch_up)
check_overlap("Stage1_UP", s1_up, "PostChemo_UP", ch_up)

cat("\n--- DOWN-REGÜLE GEN KÜMELERİ ARASI KESİŞİM (OVERLAP) ---\n")
check_overlap("Stage1_DOWN", s1_down, "Stage2_DOWN", s2_down)
check_overlap("Stage2_DOWN", s2_down, "PostChemo_DOWN", ch_down)
check_overlap("Stage1_DOWN", s1_down, "PostChemo_DOWN", ch_down)

deseq2_df <- res_Stg2_vs_Ctrl

# ==============================================================================
# ADIM 3: edgeR (QUASI-LIKELIHOOD GLM) İLE ANALİZ
# ==============================================================================
dge <- DGEList(counts = counts_matrix, group = col_data$Condition)

keep_edger <- filterByExpr(dge, group = col_data$Condition)
dge <- dge[keep_edger, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

design <- model.matrix(~ 0 + Condition, data = col_data)
colnames(design) <- levels(col_data$Condition)

dge <- estimateDisp(dge, design)
plotBCV(dge, main = "edgeR: Biyolojik Varyasyon Katsayısı (BCV)")

fit_edger <- glmQLFit(dge, design)

contrasts_edger <- makeContrasts(
  Stg1_vs_Ctrl  = Stage1 - Control,
  Stg2_vs_Ctrl  = Stage2 - Control,
  Chemo_vs_Stg2 = PostChemo - Stage2,
  levels = design
)

qlf_Stg1_vs_Ctrl  <- glmQLFTest(fit_edger, contrast = contrasts_edger[, "Stg1_vs_Ctrl"])
qlf_Stg2_vs_Ctrl  <- glmQLFTest(fit_edger, contrast = contrasts_edger[, "Stg2_vs_Ctrl"])
qlf_Chemo_vs_Stg2 <- glmQLFTest(fit_edger, contrast = contrasts_edger[, "Chemo_vs_Stg2"])

edger_res_Stg2 <- topTags(qlf_Stg2_vs_Ctrl, n = Inf)$table

summary(decideTests(qlf_Stg1_vs_Ctrl, p.value = 0.05, lfc = 1))
summary(decideTests(qlf_Stg2_vs_Ctrl, p.value = 0.05, lfc = 1))
summary(decideTests(qlf_Chemo_vs_Stg2, p.value = 0.05, lfc = 1))

# ==============================================================================
# ADIM 4: limma-voom İLE ANALİZ
# ==============================================================================
v <- voom(dge, design, plot = TRUE)

fit_limma <- lmFit(v, design)
fit_contrasts <- contrasts.fit(fit_limma, contrasts_edger)
fit_contrasts <- eBayes(fit_contrasts)

limma_res_Stg1  <- topTable(fit_contrasts, coef = "Stg1_vs_Ctrl", number = Inf, sort.by = "P")
limma_res_Stg2  <- topTable(fit_contrasts, coef = "Stg2_vs_Ctrl", number = Inf, sort.by = "P")
limma_res_Chemo <- topTable(fit_contrasts, coef = "Chemo_vs_Stg2", number = Inf, sort.by = "P")

# ==============================================================================
# ADIM 5: ÜÇ YÖNTEMİN KARŞILAŞTIRILMASI VE GÖRSELLEŞTİRME
# ==============================================================================
sig_genes_deseq2 <- rownames(subset(deseq2_df, padj < 0.05 & abs(log2FoldChange) > 1))
sig_genes_edger  <- rownames(subset(edger_res_Stg2, FDR < 0.05 & abs(logFC) > 1))
sig_genes_limma  <- rownames(subset(limma_res_Stg2, adj.P.Val < 0.05 & abs(logFC) > 1))

cat("DESeq2 anlamlı gen sayısı (Stage2 vs Ctrl):", length(sig_genes_deseq2), "\n")
cat("edgeR anlamlı gen sayısı  (Stage2 vs Ctrl):", length(sig_genes_edger), "\n")
cat("limma anlamlı gen sayısı  (Stage2 vs Ctrl):", length(sig_genes_limma), "\n")

common_genes <- Reduce(intersect, list(sig_genes_deseq2, sig_genes_edger, sig_genes_limma))
cat("Her üç yöntemde ortak bulunan gen sayısı:", length(common_genes), "\n")

limma_res_Chemo$DiffExpressed <- "Not Significant"
limma_res_Chemo$DiffExpressed[limma_res_Chemo$logFC > 1 & limma_res_Chemo$adj.P.Val < 0.05] <- "Up-regulated"
limma_res_Chemo$DiffExpressed[limma_res_Chemo$logFC < -1 & limma_res_Chemo$adj.P.Val < 0.05] <- "Down-regulated"

limma_res_Chemo$Label <- ifelse(rownames(limma_res_Chemo) %in% c(known_oncogenes, known_chemo),
                                rownames(limma_res_Chemo), NA)

ggplot(limma_res_Chemo, aes(x = logFC, y = -log10(adj.P.Val), color = DiffExpressed)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Down-regulated" = "#2b5c8f", "Not Significant" = "grey75", "Up-regulated" = "#d95f02")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text(aes(label = Label), vjust = -0.8, size = 3.2, color = "black", na.rm = TRUE) +
  labs(
    title = "Volcano Plot: FOLFOX Kemoterapisi Sonrası vs Evre II Tümör (10.000 İnsan Geni)",
    x = "log2 Kat Değişimi (log2 Fold Change)",
    y = "-log10(Düzeltilmiş P Değeri)"
  ) +
  theme_minimal()

biomarkers <- c(known_oncogenes, known_tsg, known_invasion, known_chemo)
mat_heatmap <- assay(vsd)[biomarkers, ]
mat_scaled <- t(scale(t(mat_heatmap)))

annotation_col <- data.frame(Condition = col_data$Condition)
rownames(annotation_col) <- colnames(mat_scaled)

pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  main = "Evreler ve Kemoterapi Boyunca Anahtar Biyobelirteçlerin İfadesi (Z-Score)"
)

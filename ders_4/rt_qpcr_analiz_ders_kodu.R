# ==============================================================================
# DERS 4 — RT-qPCR ANALİZİ (ÖĞRENCİ İÇİN ADIM ADIM)
# Her bölümü ayrı ayrı kopyalayıp çalıştırabilirsiniz.
# ==============================================================================

# --- 0. Paketler ve klasör ----------------------------------------------------
gerekli_paketler <- c("pcr", "rtpcr", "ggplot2", "dplyr", "tidyr", "pheatmap", "ggrepel")
yeni <- gerekli_paketler[!(gerekli_paketler %in% installed.packages()[, "Package"])]
if (length(yeni) > 0) install.packages(yeni, dependencies = TRUE)

library(pcr)
library(rtpcr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggrepel)

out_dir <- "ders_4_plot"
dir.create(out_dir, showWarnings = FALSE)
cat("Çıktı klasörü hazır:", out_dir, "\n\n")

# ==============================================================================
# BÖLÜM 1 — CONTROL vs TUMOR
# ==============================================================================

# --- 1.1 Veriyi oku ve incele -------------------------------------------------
df_tumor <- read.csv("ornek1_control_vs_tumor_qpcr.csv", stringsAsFactors = FALSE)
cat("Veri boyutları:", dim(df_tumor), "\n")
print(head(df_tumor, 3))
cat("\nGruplar:\n")
print(table(df_tumor$Group))

# --- 1.2 Ham Ct kalite kontrol grafiği ----------------------------------------
df_tumor_long <- df_tumor %>%
  pivot_longer(cols = -c(Sample_ID, Group, Replicate),
               names_to = "Gene", values_to = "Ct") %>%
  mutate(
    Gene_Type = ifelse(Gene %in% c("GAPDH", "ACTB"), "Referans Gen", "Hedef Gen"),
    Gene = factor(Gene, levels = c("GAPDH", "ACTB", "MYC", "EGFR", "VEGFA", "CCND1",
                                   "BCL2", "MMP9", "HIF1A", "TP53", "PTEN", "CDKN1A"))
  )

p_qc_ct <- ggplot(df_tumor_long, aes(x = Gene, y = Ct, fill = Group)) +
  geom_boxplot(width = 0.65, alpha = 0.88, outlier.shape = NA, color = "#1C2833") +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.65),
              size = 1.6, alpha = 0.75, color = "#1C2833") +
  scale_fill_manual(values = c("Control" = "#5D6D7E", "Tumor" = "#C0392B")) +
  labs(title = "Örnek 1: Ham Ct Değerleri (Kalite Kontrol)",
       subtitle = "Referans genler sabit kalmalı",
       x = "Genler", y = "Ct", fill = "Grup") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top")

print(p_qc_ct)
ggsave(file.path(out_dir, "01_qc_ct_boxplot.png"), p_qc_ct, width = 10, height = 6, dpi = 150)
cat("→ Grafik kaydedildi: 01_qc_ct_boxplot.png\n\n")

# --- 1.3 ΔΔCt analizi (pcr paketi) --------------------------------------------
grup_vektoru <- df_tumor$Group
ct_matris_10hedef <- df_tumor %>%
  select(GAPDH, MYC, EGFR, VEGFA, CCND1, BCL2, MMP9, TP53, PTEN, HIF1A, CDKN1A)

cat("--- pcr_analyze() sonucu ---\n")
ddct_sonuc <- pcr_analyze(
  df = ct_matris_10hedef,
  group_var = grup_vektoru,
  reference_gene = "GAPDH",
  reference_group = "Control",
  method = "delta_delta_ct",
  plot = FALSE
)
print(ddct_sonuc)

# --- 1.4 İstatistiksel test ---------------------------------------------------
cat("\n--- pcr_test() sonucu (t-test + FDR) ---\n")
test_ttest <- pcr_test(
  df = ct_matris_10hedef,
  group_var = grup_vektoru,
  reference_gene = "GAPDH",
  reference_group = "Control",
  test = "t.test"
) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    log2FC = -estimate,
    Fold_Change = 2^(-estimate),
    Regulation = ifelse(log2FC > 0, "Up-Regüle (%75)", "Down-Regüle (%25)"),
    Sig_Label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )
print(test_ttest)

# --- 1.5 Barplot (doğru SE ile) -----------------------------------------------
hedef_genler_tumor <- c("MYC", "EGFR", "VEGFA", "CCND1", "BCL2", "MMP9",
                        "TP53", "PTEN", "HIF1A", "CDKN1A")

tumor_summary <- df_tumor %>%
  select(Sample_ID, Group, GAPDH, all_of(hedef_genler_tumor)) %>%
  pivot_longer(cols = all_of(hedef_genler_tumor), names_to = "gene", values_to = "Ct") %>%
  mutate(dCt = Ct - GAPDH) %>%
  group_by(gene) %>%
  mutate(Ctrl_dCt_Mean = mean(dCt[Group == "Control"], na.rm = TRUE)) %>%
  ungroup() %>%
  filter(Group == "Tumor") %>%
  mutate(ddCt = dCt - Ctrl_dCt_Mean, log2FC = -ddCt) %>%
  group_by(gene) %>%
  summarise(
    mean_log2FC = mean(log2FC, na.rm = TRUE),
    se_log2FC = sd(log2FC, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  left_join(test_ttest %>% select(gene, Fold_Change, Regulation, Sig_Label), by = "gene") %>%
  mutate(gene = reorder(gene, mean_log2FC))

cat("\n--- Barplot için özet tablo ---\n")
print(tumor_summary)

p_tumor_bar <- ggplot(tumor_summary, aes(x = gene, y = mean_log2FC, fill = Regulation)) +
  geom_col(width = 0.68, color = "#1C2833", linewidth = 0.4, alpha = 0.92) +
  geom_errorbar(aes(ymin = mean_log2FC - se_log2FC, ymax = mean_log2FC + se_log2FC),
                width = 0.25, linewidth = 0.7, color = "#1C2833") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "#555555") +
  geom_text(aes(y = ifelse(mean_log2FC > 0,
                           mean_log2FC + se_log2FC + 0.22,
                           mean_log2FC - se_log2FC - 0.22),
                label = paste0(round(Fold_Change, 2), "x (", Sig_Label, ")")),
            size = 3.5, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = c("Up-Regüle (%75)" = "#C0392B",
                               "Down-Regüle (%25)" = "#1F618D")) +
  scale_y_continuous(limits = c(-3.2, 3.6)) +
  labs(title = "Control vs Tumor: log2(Fold Change) ± SE",
       subtitle = "Referans: GAPDH | n = 6",
       x = "Hedef Genler", y = "log2(Fold Change)", fill = "Yönelim") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        legend.position = "bottom")

print(p_tumor_bar)
ggsave(file.path(out_dir, "02_tumor_bar_plot.png"), p_tumor_bar, width = 10, height = 8, dpi = 150)
cat("→ Grafik kaydedildi: 02_tumor_bar_plot.png\n\n")

# --- 1.6 Volcano plot ---------------------------------------------------------
p_tumor_volcano <- ggplot(test_ttest, aes(x = log2FC, y = -log10(p_adj), color = Regulation)) +
  geom_point(size = 4.2, alpha = 0.9) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#555555") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "#555555") +
  geom_label_repel(aes(label = paste0(gene, " (", round(Fold_Change, 1), "x)")),
                   size = 3.4, fontface = "bold", box.padding = 0.4, show.legend = FALSE) +
  scale_color_manual(values = c("Up-Regüle (%75)" = "#C0392B",
                                "Down-Regüle (%25)" = "#1F618D")) +
  labs(title = "Volcano Plot: Tumor vs Control",
       x = "log2(Fold Change)", y = "-log10(FDR p_adj)", color = "Regülasyon") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

print(p_tumor_volcano)
ggsave(file.path(out_dir, "03_tumor_volcano.png"), p_tumor_volcano, width = 9, height = 7, dpi = 150)
cat("→ Grafik kaydedildi: 03_tumor_volcano.png\n\n")

# ==============================================================================
# BÖLÜM 2 — BALIK SİNDİRİM ENZİMLERİ
# ==============================================================================

# --- 2.1 Veriyi oku -----------------------------------------------------------
df_fish <- read.csv("ornek2_balik_sindirim_enzimleri_qpcr.csv", stringsAsFactors = FALSE)
cat("Balık veri boyutları:", dim(df_fish), "\n")
print(table(df_fish$Feed_Formulation, df_fish$Tissue))

hedef_genler_24 <- c(
  "pga1", "pga2", "tryp1", "tryp2", "chym_b", "ela1", "cpa1", "cpb1",
  "anpep", "alpi", "mgam", "si", "cel", "pnlip", "pla2g1b", "cyp7a1", "gck", "slc15a1",
  "amy2a", "chia", "lipg", "g6pc", "fhl1_est", "cath_d"
)

# --- 2.2 Faceted barplot ------------------------------------------------------
df_fish_long <- df_fish %>%
  mutate(Ref_Mean = 0.5 * (ef1a + actb)) %>%
  pivot_longer(cols = all_of(hedef_genler_24), names_to = "Gene", values_to = "Ct") %>%
  mutate(dCt = Ct - Ref_Mean)

ctrl_means <- df_fish_long %>%
  filter(Feed_Formulation == "Control_Diet") %>%
  group_by(Tissue, Gene) %>%
  summarise(Ctrl_dCt_Mean = mean(dCt), .groups = "drop")

fish_summary_df <- df_fish_long %>%
  filter(Feed_Formulation == "Formulated_Diet") %>%
  left_join(ctrl_means, by = c("Tissue", "Gene")) %>%
  mutate(ddCt = dCt - Ctrl_dCt_Mean, log2FC = -ddCt) %>%
  group_by(Tissue, Gene) %>%
  summarise(mean_log2FC = mean(log2FC),
            se_log2FC = sd(log2FC) / sqrt(n()),
            .groups = "drop") %>%
  mutate(
    Regulation = ifelse(mean_log2FC > 0, "Up-Regüle (18 Gen = %75)", "Down-Regüle (6 Gen = %25)"),
    Gene = factor(Gene, levels = hedef_genler_24),
    Tissue_TR = factor(Tissue,
                       levels = c("Stomach", "Small_Intestine", "Large_Intestine", "Liver"),
                       labels = c("Mide (Stomach)", "İnce Bağırsak (Small Intestine)",
                                  "Kalın Bağırsak (Large Intestine)", "Karaciğer (Liver)"))
  )

cat("\n--- Balık özet (ilk 8 satır) ---\n")
print(head(fish_summary_df, 8))

p_fish_facet <- ggplot(fish_summary_df, aes(x = Gene, y = mean_log2FC, fill = Regulation)) +
  geom_col(width = 0.72, color = "#1C2833", linewidth = 0.35, alpha = 0.92) +
  geom_errorbar(aes(ymin = mean_log2FC - se_log2FC, ymax = mean_log2FC + se_log2FC),
                width = 0.25, linewidth = 0.5, color = "#1C2833") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "#555555") +
  facet_wrap(~ Tissue_TR, ncol = 2) +
  scale_fill_manual(values = c("Up-Regüle (18 Gen = %75)" = "#C0392B",
                               "Down-Regüle (6 Gen = %25)" = "#1F618D")) +
  labs(title = "Balık: Yem Formülasyonuna Göre 24 Enzim Geni",
       subtitle = "Formulated_Diet vs Control_Diet | ef1a & actb",
       x = "Genler", y = "log2(Fold Change) ± SE", fill = "Yön") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "#1A365D"),
        strip.text = element_text(color = "white", face = "bold"),
        axis.text.x = element_text(angle = 50, hjust = 1, size = 8),
        legend.position = "bottom")

print(p_fish_facet)
ggsave(file.path(out_dir, "04_fish_facet_log2fc.png"), p_fish_facet, width = 14, height = 10, dpi = 150)
cat("→ Grafik kaydedildi: 04_fish_facet_log2fc.png\n\n")

# --- 2.3 Heatmap --------------------------------------------------------------
heatmap_mat <- df_fish_long %>%
  mutate(neg_dCt = -dCt,
         Grup_Etiketi = paste0(Tissue, "_", Feed_Formulation)) %>%
  group_by(Gene, Grup_Etiketi) %>%
  summarise(Mean_neg_dCt = mean(neg_dCt), .groups = "drop") %>%
  pivot_wider(names_from = Grup_Etiketi, values_from = Mean_neg_dCt) %>%
  as.data.frame()

rownames(heatmap_mat) <- heatmap_mat$Gene
heatmap_mat$Gene <- NULL
heatmap_mat <- as.matrix(heatmap_mat[hedef_genler_24, ])

gen_anotasyon <- data.frame(
  Regulasyon = factor(c(rep("Up (%75)", 18), rep("Down (%25)", 6))),
  Fonksiyon = factor(c(
    rep("Proteaz/Peptidaz", 8),
    rep("Fırçamsı Kenar/Emilim", 4),
    rep("Lipaz", 3),
    rep("Metabolizma/Taşıma", 3),
    rep("Baskılanan Enzimler", 6)
  )),
  row.names = hedef_genler_24
)

cat("--- Heatmap matrisi (ilk 5 gen, ilk 4 sütun) ---\n")
print(round(heatmap_mat[1:5, 1:4], 2))

pheatmap(
  heatmap_mat,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_row = gen_anotasyon,
  color = colorRampPalette(c("#1F618D", "white", "#C0392B"))(100),
  main = "Z-Skor İfade Haritası (-ΔCt)",
  fontsize_row = 9,
  fontsize_col = 9,
  border_color = "white",
  filename = file.path(out_dir, "05_fish_heatmap.png"),
  width = 14,
  height = 10
)
cat("→ Heatmap kaydedildi: 05_fish_heatmap.png\n\n")

cat("=== Tüm figürler kaydedildi ===\n")
cat("Klasör:", out_dir, "\n")
# Toparslan, Karabag & Bilge (2020) PLOS ONE 15(12): e0243927
# https://doi.org/10.1371/journal.pone.0243927
#
# S1 ek kodunun duzeltilmis hali. Makaledeki Fig 1-7 uretir.
# S1'deki sozdizimi hatasi, pasta eslesmesi, Kumluca/Phaselis
# yer degistirmesi ve Fig 4'un tek panelde birlestirilmemesi giderildi.
#
# Haplotip aginin dugum koordinatlari pegas tarafindan her calistirmada
# yeniden konur; makaledeki yay sekli elle oturtulmustur. Topoloji,
# daire buyuklugu, mutasyon cizgileri ve renkler makale ile aynidir.
# Fig 7 bootstrap renkleri set.seed ile sabitlenir; makale tohum
# kullanmadigi icin dugum renkleri birebir ayni olmayabilir.
#
# Calistirma: bu dosyayi S2_Appendix.fas ile ayni klasore koyun,
# sonra RStudio'da source() edin veya:
#   Rscript plosone_0243927_figures.R

fname <- "S2_Appendix.fas"
out_dir <- "figures"
dir.create(out_dir, showWarnings = FALSE)

library(ape)
library(ggplot2)
library(pegas)
library(ggtree)
library(Biostrings)

# ------------------------------------------------------------------
# 1. Oku. S2 zaten 373 bp ve hizali; msa sadece gerekirse calisir.
# ------------------------------------------------------------------
if (!file.exists(fname)) {
  stop("FASTA bulunamadi: ", fname, "\nCalisma dizinini dosyanin oldugu klasor yapin.")
}

seq_set <- readDNAStringSet(fname)
seq_chr <- toupper(gsub("[^ACGTacgtNnRYSWKMBDHV-]", "", as.character(seq_set)))
names(seq_chr) <- names(seq_set)
if (is.null(names(seq_chr)) || any(!nzchar(names(seq_chr)))) {
  stop("FASTA basliklari okunamadi. Dosya >isim ile baslamali.")
}

widths <- nchar(seq_chr)
already_aligned <- length(unique(widths)) == 1L && !any(grepl("-", seq_chr))

if (!already_aligned) {
  stop("Diziler esit uzunlukta degil. Bu lab dosyasi hizali olmali; msa kullanilmiyor.")
}
message("Diziler esit uzunlukta ve bosluksuz; hizalama atlandi (", unique(widths), " bp).")

# as.DNAbin(character) liste dondurur, rownames atanamaz. Matris olarak kur.
nbin <- as.DNAbin(do.call(rbind, strsplit(seq_chr, "")))
rownames(nbin) <- gsub(" +", "_", names(seq_chr))
rownames(nbin) <- gsub("_+", "_", rownames(nbin))

pop_of <- function(x) sub("_.*$", "", x)
pops <- pop_of(rownames(nbin))
if (any(!nzchar(pops))) stop("Populasyon adi okunamadi. Isimler Pop_no biciminde olmali.")

group_of <- function(pop) {
  ifelse(pop == "Firm", "Firm",
         ifelse(pop %in% c("Aksu", "Demre", "Kumluca"), "Greenhouse", "Nature"))
}
groups <- group_of(pops)

message("Ornek sayisi: ", nrow(nbin), "  uzunluk: ", ncol(nbin))
print(table(pops))

# Makaledeki lejant sirasi (Fig 4a)
hapcol <- c("Aksu", "Demre", "Kumluca", "Firm",
            "Bayatbadem", "Geyikbayir", "Phaselis", "Termessos")
missing_pops <- setdiff(hapcol, unique(pops))
if (length(missing_pops)) {
  warning("Beklenen populasyon yok: ", paste(missing_pops, collapse = ", "))
  hapcol <- hapcol[hapcol %in% unique(pops)]
}

# Fig 4a renkleri (makaledeki mavi / kirmizi / yesil gruplari)
pop_cols <- c(
  Aksu = "dodgerblue4", Demre = "royalblue2", Kumluca = "skyblue1",
  Firm = "red", Bayatbadem = "olivedrab4", Geyikbayir = "olivedrab3",
  Phaselis = "olivedrab1", Termessos = "darkseagreen1"
)
grp_cols <- c(Firm = "red", Greenhouse = "blue", Nature = "green")

# ------------------------------------------------------------------
# 2. NJ agaci (K80) — Fig 1, 5, 6
# ------------------------------------------------------------------
dnbin <- dist.dna(nbin, model = "K80")
tree <- nj(dnbin)

# Fig 1. Dikdortgen NJ + hizalama seridi
ggt <- ggtree(tree, aes(color = branch.length), linewidth = 0.4) +
  scale_color_continuous(high = "lightskyblue1", low = "coral4", name = "branch.length") +
  geom_tiplab(align = TRUE, size = 1.6, linesize = 0.15) +
  geom_treescale(y = -5, color = "coral4", fontsize = 3, width = 0.006) +
  theme(legend.position = "right")

nt_cols <- c(a = "rosybrown", c = "sienna1", g = "lightgoldenrod1", t = "lightskyblue1")
fig1 <- msaplot(
  ggt, nbin, offset = 0.009, width = 1, height = 0.5,
  color = nt_cols
)
ggsave(file.path(out_dir, "Fig1_msa_nj.pdf"), fig1, width = 14, height = 11)
ggsave(file.path(out_dir, "Fig1_msa_nj.png"), fig1, width = 14, height = 11, dpi = 150)

# Fig 5. Dairesel, populasyon rengi.
# Gruplar isimden kurulur; S1'deki 76:90 / 91:105 yer degistirmesi yok.
krp <- split(rownames(nbin), factor(pops, levels = hapcol))
fig5 <- ggtree(tree, layout = "circular", linewidth = 0.45) +
  xlim(-0.12, NA)
fig5 <- groupOTU(fig5, krp, "Populations") +
  aes(color = Populations) +
  geom_tiplab(size = 1.7, offset = 0.002) +
  geom_treescale(x = -0.1, color = "coral4", fontsize = 3, offset = 9) +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(linewidth = 2.5)))
ggsave(file.path(out_dir, "Fig5_nj_populations.pdf"), fig5, width = 8, height = 6.5)
ggsave(file.path(out_dir, "Fig5_nj_populations.png"), fig5, width = 8, height = 6.5, dpi = 150)

# Fig 6. Dairesel, dal uzunlugu rengi
fig6 <- ggtree(tree, layout = "circular", aes(color = branch.length), linewidth = 0.45) +
  xlim(-0.12, NA) +
  geom_tiplab(size = 1.7, offset = 0.002) +
  scale_color_continuous(high = "lightskyblue1", low = "coral4", name = "branch.length") +
  geom_treescale(x = -0.1, color = "coral4", fontsize = 3, offset = 9)
ggsave(file.path(out_dir, "Fig6_nj_distance.pdf"), fig6, width = 8, height = 6.5)
ggsave(file.path(out_dir, "Fig6_nj_distance.png"), fig6, width = 8, height = 6.5, dpi = 150)

# ------------------------------------------------------------------
# 3. Haplotip tablolari
# ------------------------------------------------------------------
seq_mat <- as.character(as.matrix(nbin))
seq_str <- apply(seq_mat, 1, paste, collapse = "")
seq_str <- toupper(seq_str)
hap_seq <- unique(seq_str)
hname <- paste0("H", seq_along(hap_seq))

var_cols <- which(apply(seq_mat, 2, function(col) length(unique(col)) > 1))
hap_index <- match(seq_str, hap_seq)
mat7 <- seq_mat[match(hap_seq, seq_str), var_cols, drop = FALSE]
rownames(mat7) <- hname
write.table(mat7, file.path(out_dir, "mat7_variable_sites.txt"),
            quote = FALSE, sep = "\t")

hfreq <- as.integer(table(factor(hap_index, levels = seq_along(hap_seq))))
ref <- mat7[1, ]
dot <- t(apply(mat7, 1, function(row) ifelse(row == ref, ".", row)))
cmstr4 <- data.frame(
  sequence = apply(dot, 1, paste, collapse = ""),
  n = hfreq,
  pct = round(100 * hfreq / sum(hfreq), 2),
  row.names = hname
)
write.table(cmstr4, file.path(out_dir, "cmstr4_haplotypes.txt"),
            quote = FALSE, sep = "\t", row.names = TRUE)

dhf <- sapply(hapcol, function(p) {
  tab <- table(factor(hap_index[pops == p], levels = seq_along(hap_seq)))
  as.integer(tab)
})
rownames(dhf) <- hname
write.table(dhf, file.path(out_dir, "dhf_pop_frequency.txt"),
            quote = FALSE, sep = "\t")

dhm <- as.matrix(dist.hamming(mat7))
write.table(dhm, file.path(out_dir, "dhm_hamming.txt"),
            quote = FALSE, sep = "\t")

# Fig 2. Isi haritasi. Altyazi "darkred-white" der; basilan sekil heat.colors'tur.
pdf(file.path(out_dir, "Fig2_heatmap.pdf"), width = 7, height = 7)
heatmap(dhm, scale = "none", col = heat.colors(100), keep.dendro = TRUE, symm = TRUE)
dev.off()
png(file.path(out_dir, "Fig2_heatmap.png"), width = 1400, height = 1400, res = 160)
heatmap(dhm, scale = "none", col = heat.colors(100), keep.dendro = TRUE, symm = TRUE)
dev.off()

# ------------------------------------------------------------------
# 4. Haplotip aglari — Fig 3 ve Fig 4
#    Pasta satirlari haplotip sirasina zorlanir (table() alfabetik bozar).
# ------------------------------------------------------------------
pie_by <- function(h, by) {
  ind.hap <- with(
    utils::stack(setNames(attr(h, "index"), rownames(h))),
    table(hap = ind, by = by[values])
  )
  ind.hap <- ind.hap[rownames(h), , drop = FALSE]
  ind.hap
}

build_net <- function(dna) {
  h <- pegas::haplotype(dna, strict = FALSE, trailingGapsAsN = TRUE)
  rownames(h) <- paste0("H", seq_len(nrow(h)))
  net <- tryCatch(
    haploNet(h, d = NULL, getProb = TRUE),
    error = function(e) haploNet(h, d = NULL, getProb = FALSE)
  )
  list(h = h, net = net)
}

net_obj <- build_net(nbin)
h <- net_obj$h
net <- net_obj$net

ind.hap <- pie_by(h, rownames(nbin))
pop.hap <- pie_by(h, pops)
pop.hap <- pop.hap[, hapcol, drop = FALSE]

ng <- nbin
rownames(ng) <- groups
netg_obj <- build_net(ng)
hg <- netg_obj$h
netg <- netg_obj$net
grp.hap <- pie_by(hg, groups)
grp.hap <- grp.hap[, c("Firm", "Greenhouse", "Nature"), drop = FALSE]

# Fig 3. Birey pastasi, rainbow
pdf(file.path(out_dir, "Fig3_hap_individuals.pdf"), width = 12, height = 7)
par(mar = c(1, 1, 1, 10))
plot(net, size = attr(net, "freq"), scale.ratio = 2, cex = 0.6,
     labels = TRUE, pie = ind.hap, show.mutation = 1, font = 2, fast = TRUE)
legend("topright", inset = c(-0.28, 0),
       colnames(ind.hap), fill = rainbow(ncol(ind.hap)),
       cex = 0.35, ncol = 2, bty = "n", xpd = NA)
dev.off()

# Fig 4. Populasyon (a) ve grup (b) ayni sayfada
pdf(file.path(out_dir, "Fig4_hap_pop_groups.pdf"), width = 12, height = 6)
layout(matrix(1:2, nrow = 1))
par(mar = c(2, 1, 1, 1))
plot(net, size = attr(net, "freq"), bg = pop_cols[colnames(pop.hap)],
     scale.ratio = 2, cex = 0.7, labels = TRUE, pie = pop.hap,
     show.mutation = 1, font = 2, fast = TRUE)
legend("topleft", hapcol, fill = pop_cols[hapcol],
       cex = 0.7, bty = "n", xpd = NA)
mtext("a", side = 1, line = 0.4, cex = 1.4, font = 2)

plot(netg, size = attr(netg, "freq"), bg = grp_cols[colnames(grp.hap)],
     scale.ratio = 2, cex = 0.7, labels = TRUE, pie = grp.hap,
     show.mutation = 1, font = 2, fast = TRUE)
legend("topleft", names(grp_cols), fill = grp_cols,
       cex = 0.8, bty = "n", xpd = NA)
mtext("b", side = 1, line = 0.4, cex = 1.4, font = 2)
dev.off()

png(file.path(out_dir, "Fig4_hap_pop_groups.png"), width = 2200, height = 1100, res = 160)
layout(matrix(1:2, nrow = 1))
par(mar = c(2, 1, 1, 1))
plot(net, size = attr(net, "freq"), bg = pop_cols[colnames(pop.hap)],
     scale.ratio = 2, cex = 0.7, labels = TRUE, pie = pop.hap,
     show.mutation = 1, font = 2, fast = TRUE)
legend("topleft", hapcol, fill = pop_cols[hapcol], cex = 0.7, bty = "n", xpd = NA)
mtext("a", side = 1, line = 0.4, cex = 1.4, font = 2)
plot(netg, size = attr(netg, "freq"), bg = grp_cols[colnames(grp.hap)],
     scale.ratio = 2, cex = 0.7, labels = TRUE, pie = grp.hap,
     show.mutation = 1, font = 2, fast = TRUE)
legend("topleft", names(grp_cols), fill = grp_cols, cex = 0.8, bty = "n", xpd = NA)
mtext("b", side = 1, line = 0.4, cex = 1.4, font = 2)
dev.off()

# ------------------------------------------------------------------
# 5. Fig 7. Haplotip NJ + bootstrap
# ------------------------------------------------------------------
set.seed(2020)
htre <- nj(dist.hamming(mat7))
bp <- boot.phylo(htre, mat7, B = 100, function(x) nj(dist.hamming(x)))
bp2 <- data.frame(node = seq_len(Nnode(htre)) + Ntip(htre), bootstrap = bp)
fig7 <- ggtree(htre, linewidth = 0.8) %<+% bp2 +
  geom_tiplab(size = 4) +
  geom_nodepoint(aes(fill = cut(bootstrap, c(0, 50, 70, 85, 100))),
                 shape = 21, size = 4) +
  theme_tree(legend.position = c(0.78, 0.22)) +
  scale_fill_manual(
    values = c("black", "red", "pink1", "white"),
    name = "Bootstrap Percentage (BP)",
    breaks = c("(85,100]", "(70,85]", "(50,70]", "(0,50]"),
    labels = expression(BP >= 85, 70 <= BP * "<85", 50 <= BP * "<70", BP < 50)
  ) +
  xlim(0, max(htre$edge.length, na.rm = TRUE) * 1.4)
ggsave(file.path(out_dir, "Fig7_haplotype_bootstrap.pdf"), fig7, width = 11, height = 6)
ggsave(file.path(out_dir, "Fig7_haplotype_bootstrap.png"), fig7, width = 11, height = 6, dpi = 150)

message("Bitti. Ciktilar: ", normalizePath(out_dir))

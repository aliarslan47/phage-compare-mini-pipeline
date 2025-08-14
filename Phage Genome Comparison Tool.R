# =========================  PHAGE COMPARISON MINI-PIPELINE  =========================
# Girdiler: T4 ve T7 için FASTA + GFF (GenBank opsiyonel)
# Çıktılar: results/ altında tablo ve görseller (barplot, circular, synteny)
# Gereken R paketleri: Biostrings, rtracklayer, seqinr, genoPlotR, circlize, tidyverse
# Dış araç GEREKMEZ (BLAST/minimap2 yok).
# ===============================================================================

suppressPackageStartupMessages({
  library(Biostrings)
  library(rtracklayer)
  library(seqinr)
  library(genoPlotR)
  library(circlize)
  library(tidyverse)
  library(ggplot2)
})

# --- AYARLAR ---
data_dir <- "C:/Users/Epigenetiks/Desktop/Phage Deneme"
out_dir  <- file.path(data_dir, "results"); if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
idA <- "Enterobacteria phage T4"
idB <- "Enterobacteria phage T7"

# --- YARDIMCI FONKSİYONLAR ---
noext     <- function(x) tools::file_path_sans_ext(basename(x))
trim_head <- function(x) sub("\\s.*$", "", x)
seq_ntlen <- function(x) nchar(as.character(x), type = "chars")          # IRanges'siz nt uzunluğu

# GC% ve uzunluk
fa_stats_one <- function(fp){
  dna <- readDNAStringSet(fp)
  L   <- sum(nchar(as.character(dna), type="chars"))
  af  <- colSums(alphabetFrequency(dna, baseOnly=TRUE))
  GCp <- 100 * (af["G"] + af["C"]) / sum(af[c("A","C","G","T")])
  tibble(id=noext(fp), total_bp=as.integer(L), GC_percent=as.numeric(GCp))
}

# GFF sayaçları (CDS / gene / tRNA)
gff_counts_one <- function(gp){
  g <- import(gp)
  tvec <- tolower(as.character(mcols(g)$type))
  has_prod <- "product" %in% tolower(colnames(mcols(g)))
  prod_col <- if (has_prod) colnames(mcols(g))[tolower(colnames(mcols(g)))=="product"] else NA_character_
  prod_txt <- if (!is.na(prod_col)) tolower(as.character(mcols(g)[[prod_col]])) else character(0)
  tibble(
    id    = noext(gp),
    CDS   = sum(tvec=="cds",  na.rm=TRUE),
    Genes = sum(tvec=="gene", na.rm=TRUE),
    tRNA  = max(sum(grepl("trna", tvec, fixed=TRUE), na.rm=TRUE),
                if (length(prod_txt)) sum(grepl("trna", prod_txt, fixed=TRUE), na.rm=TRUE) else 0L)
  )
}

# Circular için GC ve kümülatif GC-skew profili
gc_skew_df <- function(seq, win= max(1000, floor(seq_ntlen(seq)/200)),
                       step=max(200,  floor(seq_ntlen(seq)/1000))){
  L <- seq_ntlen(seq)
  if (L < 10) return(tibble(pos=numeric(0), gc=numeric(0), cum=numeric(0)))
  pos <- seq(1, max(1, L - win + 1), by=step)
  if (length(pos) == 1L) pos <- c(1, max(1, L - win + 1))
  GC  <- vapply(pos, function(p){
    seg <- subseq(seq, start=p, width=min(win, L - p + 1))
    af  <- alphabetFrequency(seg, baseOnly=TRUE)
    tot <- sum(af[c("A","C","G","T")]); if (tot==0) NA_real_ else 100*(af["G"]+af["C"])/tot
  }, numeric(1))
  skew <- vapply(pos, function(p){
    seg <- subseq(seq, start=p, width=min(win, L - p + 1))
    af  <- alphabetFrequency(seg, baseOnly=TRUE); G <- af["G"]; C <- af["C"]; s <- G + C
    if (s==0) 0 else (G - C)/s
  }, numeric(1))
  tibble(pos=pos, gc=GC, cum=cumsum(skew))
}

make_circular_plot <- function(fasta_path, gff_path, out_png){
  dna <- readDNAStringSet(fasta_path); names(dna) <- trim_head(names(dna))
  seq <- replaceAmbiguities(dna[[1]], "N")
  L   <- seq_ntlen(seq)
  
  g   <- import(gff_path)
  tvec<- tolower(as.character(mcols(g)$type))
  cds <- g[tvec=="cds"]
  cds_df <- if (length(cds)) {
    data.frame(start=as.numeric(start(cds)), end=as.numeric(end(cds)),
               strand=ifelse(as.character(strand(cds))=="+",1,-1))
  } else data.frame(start=numeric(0), end=numeric(0), strand=numeric(0))
  
  prof <- gc_skew_df(seq)
  yr_gc <- range(prof$gc, na.rm=TRUE); if (!all(is.finite(yr_gc))) yr_gc <- c(0,100)
  yr_sk <- range(prof$cum, na.rm=TRUE); if (!all(is.finite(yr_sk))) yr_sk <- c(-1,1)
  
  png(out_png, width=1400, height=1400, res=150)
  circos.clear(); circos.par(gap.after=2, start.degree=90)
  circos.initialize(factors="genome", xlim=c(0, L))
  
  circos.trackPlotRegion("genome", ylim=c(0,1), track.height=0.05,
                         panel.fun=function(x,y) circos.axis(labels.cex=0.6))
  
  if (nrow(cds_df)) {
    bed <- data.frame(chr="genome",
                      start=pmax(0, pmin(L, cds_df$start)),
                      end  =pmax(0, pmin(L, cds_df$end)),
                      strand=cds_df$strand)
    circos.genomicTrack(bed, ylim=c(0,1), track.height=0.12,
                        panel.fun=function(region, value, ...){
                          colv <- ifelse(value$strand==1, "#2C7BB6", "#D7191C")
                          circos.genomicRect(region, ybottom=0, ytop=1, col=colv, border=NA, ...)
                        })
  }
  circos.trackPlotRegion("genome", ylim=yr_gc, track.height=0.12,
                         panel.fun=function(x,y) circos.lines(prof$pos, prof$gc, lwd=1.4))
  circos.trackPlotRegion("genome", ylim=yr_sk, track.height=0.12,
                         panel.fun=function(x,y) circos.lines(prof$pos, prof$cum, lwd=1.4))
  title(paste0(noext(fasta_path), " — Circular genome (CDS, GC%, cumulative GC-skew)"))
  dev.off()
}

# AA 6-çerçeve 5-mer anchor (IRanges’siz, BLAST’sız)
GCODE <- Biostrings::getGeneticCode("11")
aa_from_frame_char <- function(seq, offset = 0){
  s <- toupper(as.character(seq)); L <- nchar(s, type="chars")
  if (L < offset + 3) return(character(0))
  n_cod <- (L - offset) %/% 3; if (n_cod <= 0) return(character(0))
  starts <- offset + 1L + (0:(n_cod - 1L)) * 3L
  codons <- substring(s, starts, starts + 2L)
  valid  <- grepl("^[ACGT]{3}$", codons)
  aa <- rep(NA_character_, n_cod); aa[valid] <- GCODE[codons[valid]]
  aa[aa=="*"] <- NA_character_; aa
}
aa_kmers <- function(aa_vec, k){
  n <- length(aa_vec); if (n < k) return(list(mers=character(0), pos=integer(0)))
  ok <- !is.na(aa_vec); cs <- cumsum(as.integer(ok))
  sums <- cs[k:n] - c(0L, cs[1:(n-k)]); idx <- which(sums == k)
  if (!length(idx)) return(list(mers=character(0), pos=integer(0)))
  mers <- vapply(idx, function(i) paste0(aa_vec[i:(i+k-1)], collapse=""), character(1))
  list(mers=mers, pos=idx)
}
map_nt_coords <- function(Lnt, frame, pos_aa, k){
  off <- as.integer(substr(frame, 2, 2))
  if (substr(frame,1,1) == "+"){
    start_nt <- off + 1L + (pos_aa - 1L) * 3L; end_nt <- start_nt + 3L*k - 1L
  } else {
    start_rc <- off + 1L + (pos_aa - 1L) * 3L; end_rc <- start_rc + 3L*k - 1L
    start_nt <- Lnt - end_rc + 1L; end_nt <- Lnt - start_rc + 1L
  }
  s <- max(1L, min(Lnt, as.integer(start_nt))); e <- max(1L, min(Lnt, as.integer(end_nt)))
  c(start=min(s,e), end=max(s,e))
}

build_anchors_aa5mer <- function(faA, faB, max_pairs=5000, keep_unique=200){
  dnaA <- readDNAStringSet(faA); names(dnaA) <- trim_head(names(dnaA))
  dnaB <- readDNAStringSet(faB); names(dnaB) <- trim_head(names(dnaB))
  seqA <- replaceAmbiguities(dnaA[[1]], "N"); LA <- seq_ntlen(seqA)
  seqB <- replaceAmbiguities(dnaB[[1]], "N"); LB <- seq_ntlen(seqB)
  
  # 6 çerçeve AA vektörleri
  rcA <- reverseComplement(seqA); rcB <- reverseComplement(seqB)
  AA_A <- list("+0"=aa_from_frame_char(seqA,0), "+1"=aa_from_frame_char(seqA,1), "+2"=aa_from_frame_char(seqA,2),
               "-0"=aa_from_frame_char(rcA,0), "-1"=aa_from_frame_char(rcA,1), "-2"=aa_from_frame_char(rcA,2))
  AA_B <- list("+0"=aa_from_frame_char(seqB,0), "+1"=aa_from_frame_char(seqB,1), "+2"=aa_from_frame_char(seqB,2),
               "-0"=aa_from_frame_char(rcB,0), "-1"=aa_from_frame_char(rcB,1), "-2"=aa_from_frame_char(rcB,2))
  
  # A için 5-mer sözlüğü
  k <- 5; cap_per_key <- 6; env <- new.env(hash=TRUE, parent=emptyenv())
  for (fr in names(AA_A)){
    km <- aa_kmers(AA_A[[fr]], k); if (!length(km$pos)) next
    spl <- split(km$pos, km$mers)
    for (key in names(spl)){
      posv <- as.integer(spl[[key]]); df_new <- data.frame(frame=fr, pos_aa=head(posv, cap_per_key))
      if (!exists(key, envir=env, inherits=FALSE)) assign(key, df_new, envir=env) else {
        old <- get(key, envir=env, inherits=FALSE); room <- cap_per_key - nrow(old[old$frame==fr,,drop=FALSE])
        if (room > 0) assign(key, rbind(old, head(df_new, room)), envir=env)
      }
    }
  }
  
  # B taraması ve anchor’lar
  anchors <- list(); idx <- 1L
  for (frB in names(AA_B)){
    km <- aa_kmers(AA_B[[frB]], k); if (!length(km$pos)) next
    for (t in seq_along(km$pos)){
      key <- km$mers[t]; if (!exists(key, envir=env, inherits=FALSE)) next
      dfA <- get(key, envir=env, inherits=FALSE); ntB <- map_nt_coords(LB, frB, km$pos[t], k)
      for (r in seq_len(nrow(dfA))){
        ntA <- map_nt_coords(LA, dfA$frame[r], dfA$pos_aa[r], k)
        anchors[[idx]] <- data.frame(
          start1=as.integer(ntA["start"]), end1=as.integer(ntA["end"]),
          start2=as.integer(ntB["start"]), end2=as.integer(ntB["end"]),
          strand_same = (substr(dfA$frame[r],1,1) == substr(frB,1,1))
        )
        idx <- idx + 1L; if (idx > max_pairs) break
      }
      if (idx > max_pairs) break
    }
    if (idx > max_pairs) break
  }
  anchors <- if (length(anchors)) do.call(rbind, anchors) else data.frame()
  if (nrow(anchors)){
    anchors$bin1 <- paste0(floor(anchors$start1/30), "_", floor(anchors$end1/30))
    anchors$bin2 <- paste0(floor(anchors$start2/30), "_", floor(anchors$end2/30))
    anchors <- anchors[ !duplicated(paste0(anchors$bin1, "|", anchors$bin2)), ]
    if (nrow(anchors) > keep_unique) anchors <- anchors[seq(1, nrow(anchors), length.out=keep_unique), ]
  }
  list(anchors=anchors, LA=LA, LB=LB)
}

make_dna_seg <- function(gff_path, Lfallback){
  g <- import(gff_path)
  tvec <- tolower(as.character(mcols(g)$type)); cds <- g[tvec=="cds"]
  if (length(cds)==0) return(as.dna_seg(data.frame(name="seg", start=1, end=Lfallback, strand=1)))
  df <- data.frame(
    name   = if (!is.null(mcols(cds)$ID)) as.character(mcols(cds)$ID) else paste0("CDS_", seq_along(cds)),
    start  = as.integer(start(cds)), end = as.integer(end(cds)),
    strand = ifelse(as.character(strand(cds))=="+",1,-1)
  )
  df <- df[order(df$start), ]; as.dna_seg(df)
}

# --- DOSYA EŞLEME ---
fa_all  <- list.files(data_dir, pattern="\\.(fa|fna|fasta)$", full.names=TRUE, ignore.case=TRUE)
gff_all <- list.files(data_dir, pattern="\\.(gff|gff3)$",     full.names=TRUE, ignore.case=TRUE)
gb_all  <- list.files(data_dir, pattern="\\.(gb|gbk|gbff)$",  full.names=TRUE, ignore.case=TRUE)

faA <- fa_all [ noext(fa_all)  == idA ][1]; faB <- fa_all [ noext(fa_all)  == idB ][1]
gffA<- gff_all[ noext(gff_all) == idA ][1]; gffB<- gff_all[ noext(gff_all) == idB ][1]
gbA <- gb_all [ noext(gb_all)  == idA ][1]; gbB <- gb_all [ noext(gb_all)  == idB ][1]
stopifnot(file.exists(faA), file.exists(faB), file.exists(gffA), file.exists(gffB))

# --- 1) TEMEL İSTATİSTİKLER ---
fa_summary <- bind_rows(fa_stats_one(faA), fa_stats_one(faB))
gff_counts <- bind_rows(gff_counts_one(gffA), gff_counts_one(gffB))
summary_tbl <- left_join(fa_summary, gff_counts, by="id")
readr::write_csv(summary_tbl, file.path(out_dir, "genome_summary_T4_T7.csv"))

# --- 2) QC BARPLOTS ---
p_len <- ggplot(summary_tbl, aes(x=reorder(id, total_bp), y=total_bp)) +
  geom_col() + coord_flip() + labs(x=NULL, y="bp", title="Genom uzunluğu")
ggsave(file.path(out_dir, "len_barplot.png"), p_len, width=6, height=4, dpi=150)

p_gc <- ggplot(summary_tbl, aes(x=reorder(id, GC_percent), y=GC_percent)) +
  geom_col() + coord_flip() + labs(x=NULL, y="GC (%)", title="GC içeriği")
ggsave(file.path(out_dir, "gc_barplot.png"), p_gc, width=6, height=4, dpi=150)

p_feat <- summary_tbl |>
  pivot_longer(cols=c(CDS, tRNA), names_to="feat", values_to="count") |>
  ggplot(aes(x=id, y=count, fill=feat)) + geom_col(position="dodge") +
  labs(x=NULL, y="Sayı", title="CDS ve tRNA sayıları") +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(out_dir, "cds_trna_barplot.png"), p_feat, width=6, height=4, dpi=150)

# --- 3) CIRCULAR HARİTALAR ---
make_circular_plot(faA, gffA, file.path(out_dir, paste0(idA, "_circular.png")))
make_circular_plot(faB, gffB, file.path(out_dir, paste0(idB, "_circular.png")))

# --- 4) ALIGNMENT-FREE K=6 COSINE (FASTA) ---
count_kmers <- function(fp, K){
  dna <- readDNAStringSet(fp); colSums(oligonucleotideFrequency(dna, width=K, step=1))
}
K <- 6
A <- count_kmers(faA, K); B <- count_kmers(faB, K)
u <- union(names(A), names(B))
A2 <- setNames(rep(0,length(u)), u); A2[names(A)] <- A
B2 <- setNames(rep(0,length(u)), u); B2[names(B)] <- B
P <- rbind(T4=A2/sum(A2), T7=B2/sum(B2)); P[is.na(P)] <- 0
L2 <- sqrt(rowSums(P^2)); L2[L2==0] <- 1e-12; Pn <- sweep(P, 1, L2, "/")
sim <- as.numeric(Pn[1,] %*% Pn[2,]); dist <- 1 - sim
write.csv(data.frame(pair="T4_vs_T7", K=K, cosine_similarity=sim, cosine_distance=dist),
          file.path(out_dir, sprintf("T4_T7_kmer%d_cosine.csv",K)), row.names=FALSE)

# --- 5) PRODUCT-BAZLI JACCARD (GFF/GB) ---
extract_products <- function(gff_path=NULL, gb_path=NULL){
  prods <- character(0)
  if (!is.null(gff_path) && file.exists(gff_path)) {
    g <- import(gff_path); cds <- g[tolower(as.character(mcols(g)$type))=="cds"]
    cols <- tolower(colnames(mcols(cds)))
    if ("product" %in% cols) {
      pcol <- colnames(mcols(cds))[cols=="product"]
      prods <- as.character(mcols(cds)[[pcol]])
    }
  }
  if ((length(prods)==0 || all(!nzchar(prods))) && !is.null(gb_path) && file.exists(gb_path)) {
    txt <- readLines(gb_path, warn=FALSE)
    hit <- grep("/product=", txt, value=TRUE)
    if (length(hit)) prods <- sub('.*?/product="?([^"]+)".*', "\\1", hit)
  }
  prods <- tolower(trimws(prods)); prods <- gsub("\\s+"," ", prods); prods <- prods[nzchar(prods)]
  bad <- c("hypothetical protein","putative protein","uncharacterized protein")
  unique(prods[!prods %in% bad])
}
setA <- extract_products(gffA, gbA); setB <- extract_products(gffB, gbB)
jac  <- function(a,b){ U <- union(a,b); if (length(U)==0) NA_real_ else length(intersect(a,b))/length(U) }
J <- jac(setA,setB)
write.csv(data.frame(pair="T4_vs_T7", jaccard=J, overlap=length(intersect(setA,setB)),
                     n_T4=length(setA), n_T7=length(setB)),
          file.path(out_dir,"T4_T7_products_jaccard.csv"), row.names=FALSE)

# --- 6) 6-FRAME AA 5-MER ANCHOR + GENOPLOTR SYNTENY ---
anc <- build_anchors_aa5mer(faA, faB, max_pairs=5000, keep_unique=200)
anchors <- anc$anchors; LA <- anc$LA; LB <- anc$LB
saveRDS(anchors, file.path(out_dir, "T4_T7_anchors_k5.rds"))

dnaA <- make_dna_seg(gffA, LA)
dnaB <- make_dna_seg(gffB, LB)
if (nrow(anchors) > 0){
  comp_df <- transform(anchors, col = ifelse(strand_same, "darkblue", "firebrick"))[
    , c("start1","end1","start2","end2","col")]
  comp <- as.comparison(comp_df)
  out_png_syn <- file.path(out_dir, "T4_vs_T7_aa5mer_genoplot.png")
  png(out_png_syn, width=1600, height=900)
  plot_gene_map(
    dna_segs = list(dnaA, dnaB),
    dna_seg_labels = c(idA, idB),
    comparisons = list(comp),
    main = "6-frame aa 5-mer anchors — T4 vs T7",
    scale = TRUE
  )
  dev.off()
}

# --- 7) ÖZET VE LİSTE ---
cat("\n== Tamam! Üretilen dosyalar ==\n")
print(list.files(out_dir, full.names=TRUE))
cat(sprintf("\nK=6 cosine similarity: %.3f (distance=%.3f)\n", sim, dist))
cat(sprintf("Product Jaccard: %.3f  (T4=%d, T7=%d, ortak=%d)\n",
            J, length(setA), length(setB), length(intersect(setA,setB))))
if (exists("out_png_syn")) cat("Synteny şekli: ", out_png_syn, "\n", sep="")
# ===============================================================================


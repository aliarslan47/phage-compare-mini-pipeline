# Phage Comparison Mini-Pipeline (T4 vs T7)

**TR ⬇️ | EN ⬇️**

---

## 🇹🇷 Türkçe

İki faj genomu için uçtan uca küçük bir karşılaştırma pipeline’ı:

- Temel istatistikler: genom uzunluğu (bp), GC(%)
- Özellik sayıları: CDS, gene, tRNA (GFF/GenBank)
- QC barplot’ları
- **Dairesel (circular) genom haritası**: CDS okları, GC%, kümülatif GC-skew
- **Alignment-free** benzerlik: k=6 k-mer cosine
- **Synteny (BLAST’sız)**: 6 çerçevede **aa 5-mer anchor** → **genoPlotR** ile çizim  
  *(BLAST/minimap2 gerekmez; tamamen R.)*

### Gereksinimler
R ≥ 4.3; paketler: `Biostrings`, `rtracklayer`, `seqinr`, `genoPlotR`, `circlize`, `tidyverse`, `ggplot2`

Kurulum:
```r
install.packages(c("tidyverse","ggplot2","seqinr","circlize","genoPlotR"))
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("Biostrings","rtracklayer"))

🇬🇧 English

An end-to-end mini pipeline to compare two phage genomes:

Basic stats: genome length (bp), GC(%)

Feature counts: CDS, gene, tRNA (from GFF/GenBank)

QC barplots

Circular genome maps: CDS arrows, GC%, cumulative GC-skew

Alignment-free similarity: k=6 k-mer cosine

Synteny (no BLAST): 6-frame aa 5-mer anchors, rendered with genoPlotR
(Runs entirely in R; no tblastx/minimap2.)

Requirements

R ≥ 4.3; packages: Biostrings, rtracklayer, seqinr, genoPlotR, circlize, tidyverse, ggplot2

Install:

install.packages(c("tidyverse","ggplot2","seqinr","circlize","genoPlotR"))
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("Biostrings","rtracklayer"))

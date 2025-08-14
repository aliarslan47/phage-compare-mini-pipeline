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

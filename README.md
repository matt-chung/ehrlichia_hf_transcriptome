# Ehrlichia sp. HF transcriptomics

# Table of Contents
<!-- MarkdownTOC autolink="true" levels="1,2,3,4" -->

- [Conduct species delineation analyses](#conduct-species-delineation-analyses)
    - [Create accession list](#create-accession-list)
  - [Compare ANI values between the Ehrlichia genomes](#compare-ani-values-between-the-ehrlichia-genomes)
    - [Calculate ANI values using OrthoANIu](#calculate-ani-values-using-orthoaniu)
    - [Process OrthoANIu output into a sequence identity matrix](#process-orthoaniu-output-into-a-sequence-identity-matrix)
    - [Plot heatmap for ANI analysis](#plot-heatmap-for-ani-analysis)
      - [Set R inputs](#set-r-inputs)
      - [Load R functions](#load-r-functions)
      - [Load R packages and view sessionInfo](#load-r-packages-and-view-sessioninfo)
      - [Construct ANI sequence identity matrix](#construct-ani-sequence-identity-matrix)
      - [Set row and column order](#set-row-and-column-order)
      - [Plot ANI heatmap](#plot-ani-heatmap)
  - [Compare dDDH values between the Ehrlichia genomes](#compare-dddh-values-between-the-ehrlichia-genomes)
    - [Calculate dDDH values using the online webservice http://ggdc.dsmz.de/](#calculate-dddh-values-using-the-online-webservice-httpggdcdsmzde)
    - [Manually construct dDDH sequence identity matrix](#manually-construct-dddh-sequence-identity-matrix)
    - [Plot heatmap for dDDH analysis](#plot-heatmap-for-dddh-analysis)
      - [Set R inputs](#set-r-inputs-1)
      - [Load R functions](#load-r-functions-1)
      - [Load R packages and view sessionInfo](#load-r-packages-and-view-sessioninfo-1)
      - [Construct dDDH sequence identity matrix](#construct-dddh-sequence-identity-matrix)
      - [Set row and column order](#set-row-and-column-order-1)
      - [Plot dDDH heatmap](#plot-dddh-heatmap)
  - [Compare CGASI values between the Ehrlichia genomes](#compare-cgasi-values-between-the-ehrlichia-genomes)
    - [Conduct core genome alignment using Mugsy](#conduct-core-genome-alignment-using-mugsy)
    - [Constructs core genome alignment fasta using only LCBs found in all genomes](#constructs-core-genome-alignment-fasta-using-only-lcbs-found-in-all-genomes)
    - [Removes all positions not present in all genomes used for core genome construction](#removes-all-positions-not-present-in-all-genomes-used-for-core-genome-construction)
    - [Plot heatmap for CGASI analysis](#plot-heatmap-for-cgasi-analysis)
      - [Set R inputs](#set-r-inputs-2)
      - [Load R functions](#load-r-functions-2)
      - [Load R packages and view sessionInfo](#load-r-packages-and-view-sessioninfo-2)
      - [Construct CGASI sequence identity matrix](#construct-cgasi-sequence-identity-matrix)
      - [Set row and column order](#set-row-and-column-order-2)
      - [Plot CGASI heatmap](#plot-cgasi-heatmap)
    - [Plot tree from core genome alignment using IQTree](#plot-tree-from-core-genome-alignment-using-iqtree)
- [Identify differentially expressed genes between HF conditions](#identify-differentially-expressed-genes-between-hf-conditions)
  - [Canine](#canine)
    - [Set R inputs](#set-r-inputs-3)
    - [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo)
    - [Create counts data frame](#create-counts-data-frame)
      - [Create TPM data frame](#create-tpm-data-frame)
      - [Set group levels](#set-group-levels)
      - [Identify differentially expressed genes longitudinally](#identify-differentially-expressed-genes-longitudinally)
  - [Tick](#tick)
    - [Set R inputs](#set-r-inputs-4)
      - [Load R functions](#load-r-functions-3)
    - [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo-1)
    - [Create counts data frame](#create-counts-data-frame-1)
    - [Create TPM data frame](#create-tpm-data-frame-1)
    - [Set group levels](#set-group-levels-1)
    - [Identify differentially expressed genes longitudinally](#identify-differentially-expressed-genes-longitudinally-1)
    - [Plot a heatmap of the up- and down-regulated genes](#plot-a-heatmap-of-the-up--and-down-regulated-genes)
      - [Create z-score log2TPM legend](#create-z-score-log2tpm-legend)
      - [Plot heatmap](#plot-heatmap)
    - [Conduct functional term enrichment analysis on up- and down-regulated gene subsets](#conduct-functional-term-enrichment-analysis-on-up--and-down-regulated-gene-subsets)

<!-- /MarkdownTOC -->

# Conduct species delineation analyses

```{bash, eval = F}
JAVA_BIN_DIR=/usr/bin

IQTREE_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/iqtree-1.6.2-Linux/bin
MOTHUR_BIN_DIR=/usr/local/packages/mothur-1.40.4
MUGSY_BIN_DIR=/local/projects-t3/EBMAL/mchung_dir/mugsy_test/mugsy_x86-64-v1r2.2
ORTHOANIU_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/OrthoANIu_v1.2
USEARCH_BIN_DIR=/usr/local/packages/usearch-9.2.64/bin
```

```{bash, eval = F}
WORKING_DIR=/local/projects-t3/EBMAL/mchung_dir/ehrlichia_hf
```

```{bash, eval = F}
mkdir -p "$WORKING_DIR"/references
mkdir -p "$WORKING_DIR"/ani
mkdir -p "$WORKING_DIR"/cgasi
mkdir -p "$WORKING_DIR"/dddh
mkdir -p "$WORKING_DIR"/plots
```

```{bash, eval = F}
wget -O "$WORKING_DIR"/references/NZ_CP007474.1.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/632/845/GCA_000632845.1_ASM63284v1/GCA_000632845.1_ASM63284v1_genomic.fna.gz
wget -O "$WORKING_DIR"/references/NC_007799.1.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/013/145/GCA_000013145.1_ASM1314v1/GCA_000013145.1_ASM1314v1_genomic.fna.gz
wget -O "$WORKING_DIR"/references/NC_023063.1.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/508/225/GCA_000508225.1_ASM50822v1/GCA_000508225.1_ASM50822v1_genomic.fna.gz
wget -O "$WORKING_DIR"/references/LANU01000000.fna.gz  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/964/755/GCA_000964755.1_ASM96475v1/GCA_000964755.1_ASM96475v1_genomic.fna.gz
wget -O "$WORKING_DIR"/references/NC_007354.1.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/012/565/GCA_000012565.1_ASM1256v1/GCA_000012565.1_ASM1256v1_genomic.fna.gz
wget -O "$WORKING_DIR"/references/NC_005295.2.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/026/005/GCA_000026005.1_ASM2600v1/GCA_000026005.1_ASM2600v1_genomic.fna.gz
wget -O "$WORKING_DIR"/references/NC_006831.1.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/050/405/GCA_000050405.1_ASM5040v1/GCA_000050405.1_ASM5040v1_genomic.fna.gz
```

```{bash, eval = F}
gunzip "$WORKING_DIR"/references/*gz
```

### Create accession list

##### Commands
```{bash, eval = F}
vim "$WORKING_DIR"/accessions.list
```

```{bash, eval = F}
NZ_CP007474.1
NC_007799.1
NC_023063.1
NZ_LANU01000001
NC_007354.1
NC_005295.2
NC_006831.1
```

## Compare ANI values between the Ehrlichia genomes

### Calculate ANI values using OrthoANIu

##### Inputs
```{bash, eval = F}
THREADS=4
FASTA_DIR="$WORKING_DIR"/references/
OUTPUT_DIR="$WORKING_DIR"/ani
```

##### Commands
```{bash, eval = F}
"$JAVA_BIN_DIR"/java -jar "$ORTHOANIU_BIN_DIR"/OAU.jar -fd "$FASTA_DIR" -fmt matrix -n "$THREADS" -o "$OUTPUT_DIR"/ani.tsv -u "$USEARCH_BIN_DIR"/usearch
```

### Process OrthoANIu output into a sequence identity matrix

##### Inputs
```{bash, eval = F}
ORTHOANIU_OUTPUT="$OUTPUT_DIR"/ani.tsv
```

##### Commands
```{bash, eval = F}
grep -B 99999 "# OrthoANIu results as matrix" "$ORTHOANIU_OUTPUT" | tail -n+3 | head -n-4 | cut -f2 | sed "s/[.]fna//g" > "$(dirname "$ORTHOANIU_OUTPUT")"/temp1
grep -A 99999 "# OrthoANIu results as matrix" "$ORTHOANIU_OUTPUT" | tail -n+3 | head -n-1 | cut -f2- > "$(dirname "$ORTHOANIU_OUTPUT")"/temp2
paste "$(dirname "$ORTHOANIU_OUTPUT")"/temp1 "$(dirname "$ORTHOANIU_OUTPUT")"/temp2 > "$ORTHOANIU_OUTPUT"
rm "$(dirname "$ORTHOANIU_OUTPUT")"/temp1
rm "$(dirname "$ORTHOANIU_OUTPUT")"/temp2
```

### Plot heatmap for ANI analysis

#### Set R inputs
```{R}
WORKING.DIR <- "Z:/EBMAL/mchung_dir/ehrlichia_hf/"
ORTHOANIU_OUTPUT.PATH <- "Z:/EBMAL/mchung_dir/ehrlichia_hf/ani/ani.tsv"
```

#### Load R functions
```{R}
formalize_sample_names <- function(x){
  x <- gsub("NZ_CP007474.1", "Ehrlichia sp. HF", x)
  x <- gsub("NC_007799.1", "E. chaffeensis Arkansas", x)
  x <- gsub("NC_023063.1", "E. muris AS145", x)
  x <- gsub("LANU01000000", "E. muris subsp. eauclairensis Wisconsin", x)
  x <- gsub("NC_007354.1", "E. canis Jake", x)
  x <- gsub("NC_005295.2", "E. ruminantium Welgevonden", x)
  x <- gsub("NC_006831.1", "E. ruminantium Gardel", x)
}
```

#### Load R packages and view sessionInfo
```{R}
require(Biostrings)
require(ggplot2)
require(pvclust)
require(reshape)

sessionInfo()
```

```{R, eval = F}
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18362)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] reshape_0.8.8       Biostrings_2.50.2   XVector_0.22.0      IRanges_2.16.0      S4Vectors_0.20.1    BiocGenerics_0.28.0
 [7] RCurl_1.95-4.12     bitops_1.0-6        pvclust_2.0-0       gridExtra_2.3       ggplot2_3.2.0       gplots_3.0.1.1     
[13] edgeR_3.24.3        limma_3.38.3        dendextend_1.12.0  

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.2         plyr_1.8.4         pillar_1.4.2       compiler_3.5.0     zlibbioc_1.28.0    viridis_0.5.1      tools_3.5.0       
 [8] digest_0.6.20      tibble_2.1.3       gtable_0.3.0       viridisLite_0.3.0  lattice_0.20-35    pkgconfig_2.0.2    rlang_0.4.0       
[15] rstudioapi_0.10    yaml_2.2.0         xfun_0.8           withr_2.1.2        dplyr_0.8.3        knitr_1.23         gtools_3.8.1      
[22] caTools_1.17.1.2   locfit_1.5-9.1     grid_3.5.0         tidyselect_0.2.5   glue_1.3.1         R6_2.4.0           gdata_2.18.0      
[29] purrr_0.3.2        magrittr_1.5       splines_3.5.0      scales_1.0.0       assertthat_0.2.1   colorspace_1.4-1   labeling_0.3      
[36] KernSmooth_2.23-15 lazyeval_0.2.2     munsell_0.5.0      crayon_1.3.4
```

#### Construct ANI sequence identity matrix

```{R}
sim <- read.delim(ORTHOANIU_OUTPUT.PATH, header = F, row.names = 1)
colnames(sim) <- rownames(sim)
```

#### Set row and column order

```{R}
rownames(sim) <- formalize_sample_names(rownames(sim))
colnames(sim) <- formalize_sample_names(colnames(sim))

levels <- pvclust(sim, nboot=5)
levels <- rownames(sim)[levels$hclust$order]
```

#### Plot ANI heatmap

```{R, fig.height=5,fig.width=6}
sim <- sim[match(levels,rownames(sim)),match(levels,colnames(sim))]
sim[lower.tri(sim)] <- NA

plot.df <- as.data.frame(cbind(rownames(sim),
                               sim))
names(plot.df)[1] <- "temp"
plot.df <- melt(plot.df,id.vars="temp", na.rm=T)
colnames(plot.df) <- c("Var1", "Var2", "value")
plot.df$Var1 <- factor(plot.df$Var1, levels=levels)
plot.df$Var2 <- factor(plot.df$Var2, levels=rev(levels))

ani.hm <- ggplot(data = plot.df, aes(Var1, Var2, fill = value))+
  geom_tile(color = "black")+
  geom_text(aes(label = ifelse(value>=0, round(value, 1), "")), size = 3)+
  scale_fill_gradientn(colors = colorRampPalette(c("blue","darkturquoise","darkgreen","green","yellow","red"))(11),
                       limits=c(75,100))+
  theme_minimal()+
  labs(x="",y="",fill="ANI")+
  #guides(fill = F)+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1,size=8),
        axis.text.y = element_text(size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pdf(paste0(WORKING.DIR,"/plots/ani.pdf"),
    height=5,
    width=6)
print(ani.hm)
dev.off()

png(paste0(WORKING.DIR,"/plots/ani.png"),
    height=5,
    width=6,
    units = "in",res=300)
print(ani.hm)
dev.off()

print(ani.hm)
```

![image](/images/ani.png)

## Compare dDDH values between the Ehrlichia genomes

### Calculate dDDH values using the online webservice http://ggdc.dsmz.de/

##### Inputs
```{bash, eval = F}
FASTA_DIR="$WORKING_DIR"/references/
OUTPUT_DIR="$WORKING_DIR"/dddh
```

##### Commands
```{bash, eval = F}
for FASTA1 in $(find "$FASTA_DIR" -name "*.fna")
do
  for FASTA2 in $(find "$FASTA_DIR" -name "*.fna")
  do
    "$JAVA_BIN_DIR"/java -jar "$ORTHOANITOOL_BIN_DIR"/OAT_cmd.jar -blastplus_dir "$NCBIBLAST_BIN_DIR" -method ggdc -fasta1 "$FASTA1" -fasta2 "$FASTA2" > "$OUTPUT_DIR"/"$(basename "$FASTA1" | sed "s/[.]fna//g")"_"$(basename "$FASTA2" | sed "s/[.]fna//g")".ggdc
  done
done
```

### Manually construct dDDH sequence identity matrix

### Plot heatmap for dDDH analysis

#### Set R inputs
```{R}
WORKING.DIR <- "Z:/EBMAL/mchung_dir/ehrlichia_hf/"
DDDH_OUTPUT.PATH <- "Z:/EBMAL/mchung_dir/ehrlichia_hf/dddh/dddh.tsv"
```

#### Load R functions
```{R}
formalize_sample_names <- function(x){
  x <- gsub("NZ_CP007474.1", "Ehrlichia sp. HF", x)
  x <- gsub("NC_007799.1", "E. chaffeensis Arkansas", x)
  x <- gsub("NC_023063.1", "E. muris AS145", x)
  x <- gsub("LANU01000000", "E. muris subsp. eauclairensis Wisconsin", x)
  x <- gsub("NC_007354.1", "E. canis Jake", x)
  x <- gsub("NC_005295.2", "E. ruminantium Welgevonden", x)
  x <- gsub("NC_006831.1", "E. ruminantium Gardel", x)
}
```

#### Load R packages and view sessionInfo

```{R}
require(Biostrings)
require(ggplot2)
require(pvclust)
require(reshape)

sessionInfo()
```

```{R, eval = F}
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18362)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] reshape_0.8.8       Biostrings_2.50.2   XVector_0.22.0      IRanges_2.16.0      S4Vectors_0.20.1    BiocGenerics_0.28.0
 [7] RCurl_1.95-4.12     bitops_1.0-6        pvclust_2.0-0       gridExtra_2.3       ggplot2_3.2.0       gplots_3.0.1.1     
[13] edgeR_3.24.3        limma_3.38.3        dendextend_1.12.0  

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.2         plyr_1.8.4         pillar_1.4.2       compiler_3.5.0     zlibbioc_1.28.0    viridis_0.5.1      tools_3.5.0       
 [8] digest_0.6.20      tibble_2.1.3       gtable_0.3.0       viridisLite_0.3.0  lattice_0.20-35    pkgconfig_2.0.2    rlang_0.4.0       
[15] rstudioapi_0.10    yaml_2.2.0         xfun_0.8           withr_2.1.2        dplyr_0.8.3        knitr_1.23         gtools_3.8.1      
[22] caTools_1.17.1.2   locfit_1.5-9.1     grid_3.5.0         tidyselect_0.2.5   glue_1.3.1         R6_2.4.0           gdata_2.18.0      
[29] purrr_0.3.2        magrittr_1.5       splines_3.5.0      scales_1.0.0       assertthat_0.2.1   colorspace_1.4-1   labeling_0.3      
[36] KernSmooth_2.23-15 lazyeval_0.2.2     munsell_0.5.0      crayon_1.3.4
```

#### Construct dDDH sequence identity matrix

```{R}
sim <- read.delim(DDDH_OUTPUT.PATH, header = F, row.names = 1)
colnames(sim) <- rownames(sim)
```

#### Set row and column order

```{R}
rownames(sim) <- formalize_sample_names(rownames(sim))
colnames(sim) <- formalize_sample_names(colnames(sim))

# levels <- pvclust(sim, nboot=5)
# levels <- rownames(sim)[levels$hclust$order]
levels <- c("E. ruminantium Welgevonden","E. ruminantium Gardel","E. canis Jake","E. chaffeensis Arkansas","Ehrlichia sp. HF", "E. muris AS145",    "E. muris subsp. eauclairensis Wisconsin")
```

#### Plot dDDH heatmap

```{R, fig.height=5,fig.width=6}
sim <- sim[match(levels,rownames(sim)),match(levels,colnames(sim))]
sim[lower.tri(sim)] <- NA

plot.df <- as.data.frame(cbind(rownames(sim),
                               sim))
names(plot.df)[1] <- "temp"
plot.df <- melt(plot.df,id.vars="temp", na.rm=T)
colnames(plot.df) <- c("Var1", "Var2", "value")
plot.df$Var1 <- factor(plot.df$Var1, levels=levels)
plot.df$Var2 <- factor(plot.df$Var2, levels=rev(levels))

dddh.hm <- ggplot(data = plot.df, aes(Var1, Var2, fill = value))+
  geom_tile(color = "black")+
  geom_text(aes(label = ifelse(value>=0, round(value, 1), "")), size = 3)+
  scale_fill_gradientn(colors = colorRampPalette(c("blue","darkturquoise","darkgreen","green","yellow","red"))(11),
                       limits=c(20,100))+
  theme_minimal()+
  labs(x="",y="",fill="dDDH")+
  #guides(fill = F)+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1,size=8),
        axis.text.y = element_text(size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pdf(paste0(WORKING.DIR,"/plots/dddh.pdf"),
    height=5,
    width=6)
print(dddh.hm)
dev.off()

png(paste0(WORKING.DIR,"/plots/dddh.png"),
    height=5,
    width=6,
    units = "in",res=300)
print(dddh.hm)
dev.off()

print(dddh.hm)
```

![image](/images/dddh.png)

## Compare CGASI values between the Ehrlichia genomes

### Conduct core genome alignment using Mugsy

##### Inputs
```{bash, eval = F}
FASTA_DIR="$WORKING_DIR"/references/
OUTPUT_DIR="$WORKING_DIR"/cgasi/
```

##### Commands
```{bash, eval = F}
source "$MUGSY_BIN_DIR"/mugsyenv.sh
"$MUGSY_BIN_DIR"/mugsy --prefix mugsy --directory "$OUTPUT_DIR" $(find "$FASTA_DIR" -name "*[.]fna")
```

### Constructs core genome alignment fasta using only LCBs found in all genomes

##### Inputs
```{bash, eval = F}
ACCESSIONS_LIST="$WORKING_DIR"/accessions.list
MUGSY_OUTPUT_DIR="$WORKING_DIR"/cgasi/
```

##### Commands
```{bash, eval = F}
mkdir "$MUGSY_OUTPUT_DIR"/processed/
rm "$MUGSY_OUTPUT_DIR"/processed/*
while read LINE
do
  echo ">"$LINE"" > "$MUGSY_OUTPUT_DIR"/processed/"$(echo "$LINE" | sed "s/.[.1-9]$//g")".fna
done < "$ACCESSIONS_LIST"

while read LINE
do
  MARKER=$(echo "$LINE" | awk -F " " '{print $1}')
  if [ "$MARKER" = "a" ]
  then
    MULT=$(echo "$LINE" | awk -F " " '{print $4}')
  else
    if [ "$MULT" = "mult="$(wc -l "$ACCESSIONS_LIST" | awk '{print $1}')"" ]
    then
      SPECIES=$(echo "$LINE" | awk -F " " '{print $2}' | sed "s/[.].*//g")
      echo "$LINE" | awk -F " " '{print $7}' >> "$MUGSY_OUTPUT_DIR"/processed/"$SPECIES".fna
    fi
  fi
done < "$MUGSY_OUTPUT_DIR"/mugsy.maf
```

### Removes all positions not present in all genomes used for core genome construction

##### Inputs
```{bash, eval = F}
MUGSY_OUTPUT_DIR="$WORKING_DIR"/cgasi/
```

##### Commands
```{bash, eval = F}
cat "$MUGSY_OUTPUT_DIR"/processed/* > "$MUGSY_OUTPUT_DIR"/concat.fasta
"$MOTHUR_BIN_DIR"/mothur #filter.seqs(fasta="$MUGSY_OUTPUT_DIR"/concat.fasta, vertical=F, trump=-)
"$MOTHUR_BIN_DIR"/mothur "#filter.seqs(fasta="$MUGSY_OUTPUT_DIR"/concat.filter.fasta, vertical=F, trump=.)"
mv "$MUGSY_OUTPUT_DIR"/concat.filter.filter.fasta "$MUGSY_OUTPUT_DIR"/concat_alignment.fna
```

### Plot heatmap for CGASI analysis

Core genome alignment is 530,738 bp long.

#### Set R inputs
```{R}
WORKING.DIR <- "Z:/EBMAL/mchung_dir/ehrlichia_hf/"
CONCAT_ALIGNMENT.PATH <- "Z:/EBMAL/mchung_dir/ehrlichia_hf/cgasi/concat_alignment.fna"
```

#### Load R functions
```{R}
formalize_sample_names <- function(x){
  x <- gsub("NZ_CP007474.1", "Ehrlichia sp. HF", x)
  x <- gsub("NC_007799.1", "E. chaffeensis Arkansas", x)
  x <- gsub("NC_023063.1", "E. muris AS145", x)
  x <- gsub("LANU01000000", "E. muris subsp. eauclairensis Wisconsin", x)
  x <- gsub("NC_007354.1", "E. canis Jake", x)
  x <- gsub("NC_005295.2", "E. ruminantium Welgevonden", x)
  x <- gsub("NC_006831.1", "E. ruminantium Gardel", x)
}
```

#### Load R packages and view sessionInfo
```{R}
require(Biostrings)
require(ggplot2)
require(pvclust)
require(reshape)

sessionInfo()
```

```{R, eval = F}
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18362)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] reshape_0.8.8       Biostrings_2.50.2   XVector_0.22.0      IRanges_2.16.0      S4Vectors_0.20.1    BiocGenerics_0.28.0
 [7] RCurl_1.95-4.12     bitops_1.0-6        pvclust_2.0-0       gridExtra_2.3       ggplot2_3.2.0       gplots_3.0.1.1     
[13] edgeR_3.24.3        limma_3.38.3        dendextend_1.12.0  

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.2         plyr_1.8.4         pillar_1.4.2       compiler_3.5.0     zlibbioc_1.28.0    viridis_0.5.1      tools_3.5.0       
 [8] digest_0.6.20      tibble_2.1.3       gtable_0.3.0       viridisLite_0.3.0  lattice_0.20-35    pkgconfig_2.0.2    rlang_0.4.0       
[15] rstudioapi_0.10    yaml_2.2.0         xfun_0.8           withr_2.1.2        dplyr_0.8.3        knitr_1.23         gtools_3.8.1      
[22] caTools_1.17.1.2   locfit_1.5-9.1     grid_3.5.0         tidyselect_0.2.5   glue_1.3.1         R6_2.4.0           gdata_2.18.0      
[29] purrr_0.3.2        magrittr_1.5       splines_3.5.0      scales_1.0.0       assertthat_0.2.1   colorspace_1.4-1   labeling_0.3      
[36] KernSmooth_2.23-15 lazyeval_0.2.2     munsell_0.5.0      crayon_1.3.4
```

#### Construct CGASI sequence identity matrix

```{R}
concat.alignment <- readDNAMultipleAlignment(CONCAT_ALIGNMENT.PATH)

sim <- as.data.frame(matrix(nrow = length(concat.alignment@unmasked),
                            ncol = length(concat.alignment@unmasked)))
rownames(sim) <- names(concat.alignment@unmasked)
colnames(sim) <- names(concat.alignment@unmasked)

for(i in 1:nrow(sim)){
  j <- 1
  while(j <= i){
    palign <- PairwiseAlignments(c(concat.alignment@unmasked[i],concat.alignment@unmasked[j]))
    sim[i,j] <- pid(palign,type="PID4")
    j <- j+1
  }
}

for(i in 1:nrow(sim)){
  sim[i,] <- t(sim[,i])
}
```

#### Set row and column order

```{R}
rownames(sim) <- formalize_sample_names(rownames(sim))
colnames(sim) <- formalize_sample_names(colnames(sim))

levels <- pvclust(sim, nboot=5)
levels <- rownames(sim)[levels$hclust$order]
```

#### Plot CGASI heatmap

```{R, fig.height=5,fig.width=6}
sim <- sim[match(levels,rownames(sim)),match(levels,colnames(sim))]
sim[lower.tri(sim)] <- NA

plot.df <- as.data.frame(cbind(rownames(sim),
                               sim))
names(plot.df)[1] <- "temp"
plot.df <- melt(plot.df,id.vars="temp", na.rm=T)
colnames(plot.df) <- c("Var1", "Var2", "value")
plot.df$Var1 <- factor(plot.df$Var1, levels=levels)
plot.df$Var2 <- factor(plot.df$Var2, levels=rev(levels))

cgasi.hm <- ggplot(data = plot.df, aes(Var1, Var2, fill = value))+
  geom_tile(color = "black")+
  geom_text(aes(label = ifelse(value>=0, round(value, 1), "")), size = 3)+
  scale_fill_gradientn(colors = colorRampPalette(c("blue","darkturquoise","darkgreen","green","yellow","red"))(11),
                       limits=c(80,100))+
  theme_minimal()+
  labs(x="",y="",fill="CGASI")+
  #guides(fill = F)+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1,size=8),
        axis.text.y = element_text(size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pdf(paste0(WORKING.DIR,"/plots/cgasi.pdf"),
    height=5,
    width=6)
print(cgasi.hm)
dev.off()

png(paste0(WORKING.DIR,"/plots/cgasi.png"),
    height=5,
    width=6,
    units = "in",res=300)
print(cgasi.hm)
dev.off()

print(cgasi.hm)
```

![image](/images/cgasi.png)

### Plot tree from core genome alignment using IQTree

##### Inputs
```{bash, eval = F}
THREADS=4
MSA_FNA="$WORKING_DIR"/cgasi//concat_alignment.fna
```

##### Commands
```{bash, eval = F}
"$IQTREE_BIN_DIR"/iqtree -s "$MSA_FNA" -nt "$THREADS" -bb 1000 -redo
```

# Identify differentially expressed genes between HF conditions

Analysis was not run on Ehrlichia, because samples do not have enough reads

```{R, eval = F}
            PECHA_HF_totalRNA1             PECHA_HF_totalRNA2            PECHA_HF_totalRNA3
                       1683.63                        2000.25                        931.12
       PECHA_ISE6_HF_totalRNA1        PECHA_ISE6_HF_totalRNA2       PECHA_ISE6_HF_totalRNA3    
                       5678.74                         794.11                       4446.62
```

## Canine

### Set R inputs
```{R}
WORKING.DIR <- "Z:/EBMAL/mchung_dir/PECHA/"
SALMON_OUTPUT.DIR <- "Z:/EBMAL/mchung_dir/PECHA/salmon"
GENEINFO.PATH <- "Z:/EBMAL/mchung_dir/PECHA/references/canine.cds.fna.interproscan.geneinfo.tsv"
GROUPS.PATH <- "Z:/EBMAL/mchung_dir/PECHA/pecha_groups.tsv"
```

### Load packages and view sessionInfo

```{R}
library(edgeR)

sessionInfo()
```

```{R, eval = F}
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18362)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] edgeR_3.24.3 limma_3.38.3

loaded via a namespace (and not attached):
[1] compiler_3.5.0  tools_3.5.0     Rcpp_1.0.2      grid_3.5.0      locfit_1.5-9.1  knitr_1.23      xfun_0.8        lattice_0.20-35
```

### Create counts data frame
```{R}
groups <- read.delim(GROUPS.PATH, header = F)
groups <- groups[c(grep("_DH82_",groups[,2]),grep("_HF_",groups[,2])),]
groups <- groups[intersect(grep("mRNA", groups[,2]),grep("ISE6", groups[,2], invert = T)),]

rownames <- as.character(read.delim(paste0(SALMON_OUTPUT.DIR, "/", groups[1,1],"/quant.sf"))[grep("ENSCAFT", read.delim(paste0(SALMON_OUTPUT.DIR, "/", groups[1,1],"/quant.sf"))[,1]),1])
colnames <- unique(groups[,2])

counts <- as.data.frame(matrix(0,
                               nrow = length(rownames),
                               ncol = length(colnames)))
rownames(counts) <- rownames
colnames(counts) <- colnames

for(i in 1:ncol(counts)){
  srr.vector <- groups[groups[,2] == colnames(counts[i]),1]
  for(j in 1:length(srr.vector)){
    counts.subset <- read.delim(paste0(SALMON_OUTPUT.DIR, "/",srr.vector[j],"/quant.sf"))
    counts[,i] <- counts[,i] + counts.subset[match(rownames(counts),counts.subset[,1]),5]
  }
}

colSums(counts)
```

```{R, eval = F}
PECHA_DH82_mRNA1 PECHA_DH82_mRNA2 PECHA_DH82_mRNA3
          612736          1219322          2046861
  PECHA_HF_mRNA1   PECHA_HF_mRNA2   PECHA_HF_mRNA3
         1748842          1136319          1523891
```

#### Create TPM data frame
```{R}
genelength <- read.delim(paste0(SALMON_OUTPUT.DIR, "/", groups[1,1],"/quant.sf"))
genelength <- genelength[match(rownames(counts),genelength[,1]),2]

tpm <- counts
for(i in 1:ncol(tpm)){
  tpm[,i] <- tpm[,i]/genelength
  tpm[,i] <- tpm[,i]/(sum(tpm[,i])/1000000)
}

dim(tpm)
```

```{R, eval = F}
[1] 45094     6
```

#### Set group levels
```{R}
groups <- unique(groups[2:5])
groups[,1] <- factor(groups[,1], levels = groups[,1])
groups[,2] <- factor(groups[,2], levels = unique(groups[,2]))
groups[,3] <- factor(groups[,3], levels = unique(groups[,3]))
groups[,4] <- factor(groups[,4], levels = unique(groups[,4]))
```

#### Identify differentially expressed genes longitudinally

edgeR and DESeq2 are both run with a FDR cutoff of <0.05 and a minimum CPM cutoff of 5 reads in the lowest sequenced sample in the data set.  

```{R}
FDRcutoff <- 0.05
cpm.cutoff <- 5/min(colSums(counts)) * 1000000

y <- DGEList(counts = counts, group = groups[,4])
y <- calcNormFactors(y)
keep <- rowSums(cpm(y) >= cpm.cutoff) >= min(table(groups[,4]))
keep.df <- as.data.frame(table(keep))
print(paste0(keep.df[keep.df[,1] == F,2]," genes excluded with CPM cutoff"))

y <- y[keep, , keep.lib.sizes = F]
design <- model.matrix(~groups[,4])
y <- estimateDisp(y , design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2:ncol(fit))

FDR <- as.data.frame(p.adjust(qlf$table$PValue, method="BH"))
rownames(FDR) <- rownames(qlf$table)
FDR.genes <- rownames(FDR[FDR[,1] < FDRcutoff, , drop = F])
print(paste0(length(FDR.genes)," DE genes identified using edgeR pairwise" ))

qlf$table <- as.data.frame(cbind(qlf$table,
                                 FDR))
colnames(qlf$table)[5] <- "FDR"
write.table(qlf$table[qlf$table$FDR < 0.05,],
            paste0(WORKING.DIR,"/canine_de_logFC.tsv"),
            col.names = T,
            row.names = T,
            quote = F,
            sep = "\t")
```

```{R, eval = F}
[1] "34087 genes excluded with CPM cutoff"
[1] "2 DE genes identified using edgeR pairwise
```

## Tick

### Set R inputs
```{R}
WORKING.DIR <- "Z:/EBMAL/mchung_dir/PECHA/"
SALMON_OUTPUT.DIR <- "Z:/EBMAL/mchung_dir/PECHA/salmon"
GENEINFO.PATH <- "Z:/EBMAL/mchung_dir/PECHA/references/tick.cds.fna.interproscan.geneinfo.tsv"
GROUPS.PATH <- "Z:/EBMAL/mchung_dir/PECHA/pecha_groups.tsv"
```

#### Load R functions

```{R}
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
}
functionaltermenrichment <- function(genes, geneinfo){
  for(i in 1:ncol(geneinfo)){geneinfo[,i] <- as.character(geneinfo[,i])}
  geneinfo$interpro_description[which(is.na(geneinfo$interpro_description))] <- "No InterPro entry"
  geneinfo$go_biologicalprocess[which(is.na(geneinfo$go_biologicalprocess))] <- "No GO terms for biological process"
  geneinfo$go_cellularcomponent[which(is.na(geneinfo$go_cellularcomponent))] <- "No GO terms for cellular component"
  geneinfo$go_molecularfunction[which(is.na(geneinfo$go_molecularfunction))] <- "No GO terms for molecular function"
  
  functionalterms.list <- list(ipr=as.data.frame(table(unlist(strsplit(paste(geneinfo$interpro_description, collapse = "|"),  split = "[|]")))),
                               gobio=as.data.frame(table(unlist(strsplit(paste(geneinfo$go_biologicalprocess, collapse = "|"),  split = "[|]")))),
                               gocell=as.data.frame(table(unlist(strsplit(paste(geneinfo$go_cellularcomponent, collapse = "|"),  split = "[|]")))),
                               gomol=as.data.frame(table(unlist(strsplit(paste(geneinfo$go_molecularfunction, collapse = "|"),  split = "[|]")))))
  
  geneinfo.subset <- geneinfo[geneinfo$gene %in% genes,]
  term <- c()
  clusteroccurences <- c()
  genomeoccurences <- c()
  pvalue <- c()
  correctedpvalue <- c()
  oddsratio <- c()

  functionalterms.list.subset <- list(ipr=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$interpro_description, collapse = "|"),  split = "[|]")))),
                                      gobio=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$go_biologicalprocess, collapse = "|"),  split = "[|]")))),
                                      gocell=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$go_cellularcomponent, collapse = "|"),  split = "[|]")))),
                                      gomol=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$go_molecularfunction, collapse = "|"),  split = "[|]")))))
  
  for(i in 1:length(functionalterms.list)){
    for(j in 1:nrow(functionalterms.list[[i]])){
      freq.all <- functionalterms.list[[i]][j,2]
      freq.subset <- ifelse(functionalterms.list[[i]][j,1] %in% functionalterms.list.subset[[i]][,1],
                            functionalterms.list.subset[[i]][functionalterms.list.subset[[i]][,1] == as.character(functionalterms.list[[i]][j,1]),2],
                            0)
      genes.all <- nrow(geneinfo)
      genes.subset <- nrow(geneinfo.subset)

      fisherexact.matrix <- matrix(c(freq.subset, freq.all - freq.subset,
                                     genes.subset - freq.subset, genes.all - genes.subset - freq.all + freq.subset),
                                   nrow = 2,
                                   ncol = 2)
      fisher.test <- fisher.test(fisherexact.matrix)
      
      term[length(term) + 1] <- as.character(functionalterms.list[[i]][j,1])
      clusteroccurences[length(clusteroccurences) + 1] <- as.numeric(as.character(freq.subset))
      genomeoccurences[length(genomeoccurences) + 1] <- as.numeric(as.character(freq.all))
      pvalue[length(pvalue) + 1] <- as.numeric(as.character(fisher.test$p.value))
      correctedpvalue[length(correctedpvalue) + 1] <- p.adjust(as.numeric(as.character(fisher.test$p.value)), method = "fdr", n = nrow(functionalterms.list[[i]]))
      oddsratio[length(oddsratio) + 1] <- as.numeric(as.character(fisher.test$estimate))
    }
  }
  
  terms.df <- as.data.frame(cbind(term,
                                  clusteroccurences,
                                  genomeoccurences,
                                  pvalue,
                                  correctedpvalue,
                                  oddsratio))
  terms.df <- terms.df[order(as.numeric(as.character(terms.df$pvalue))),]
  return(terms.df)
}
```

### Load packages and view sessionInfo

```{R}
library(dendextend)
library(edgeR)
library(gplots)
library(ggplot2)
library(gridExtra)
library(pvclust)

sessionInfo()
```

```{R, eval = F}
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18362)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] pvclust_2.0-0     gridExtra_2.3     ggplot2_3.2.0     gplots_3.0.1.1    edgeR_3.24.3      limma_3.38.3      dendextend_1.12.0

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.2         pillar_1.4.2       compiler_3.5.0     viridis_0.5.1      bitops_1.0-6       tools_3.5.0        tibble_2.1.3      
 [8] gtable_0.3.0       viridisLite_0.3.0  lattice_0.20-35    pkgconfig_2.0.2    rlang_0.4.0        rstudioapi_0.10    xfun_0.8          
[15] withr_2.1.2        dplyr_0.8.3        knitr_1.23         gtools_3.8.1       caTools_1.17.1.2   locfit_1.5-9.1     grid_3.5.0        
[22] tidyselect_0.2.5   glue_1.3.1         R6_2.4.0           gdata_2.18.0       purrr_0.3.2        magrittr_1.5       scales_1.0.0      
[29] assertthat_0.2.1   colorspace_1.4-1   KernSmooth_2.23-15 lazyeval_0.2.2     munsell_0.5.0      crayon_1.3.4      
```

### Create counts data frame
```{R}
groups <- read.delim(GROUPS.PATH, header = F)
groups <- groups[c(grep("_ISE6_mRNA",groups[,2]),grep("_HF_",groups[,2])),]
groups <- groups[intersect(grep("mRNA", groups[,2]),grep("ISE6", groups[,2])),]

rownames <- as.character(read.delim(paste0(SALMON_OUTPUT.DIR, "/", groups[1,1],"/quant.sf"))[grep("ISCW", read.delim(paste0(SALMON_OUTPUT.DIR, "/", groups[1,1],"/quant.sf"))[,1]),1])
colnames <- unique(groups[,2])

counts <- as.data.frame(matrix(0,
                               nrow = length(rownames),
                               ncol = length(colnames)))
rownames(counts) <- rownames
colnames(counts) <- colnames

for(i in 1:ncol(counts)){
  srr.vector <- groups[groups[,2] == colnames(counts[i]),1]
  for(j in 1:length(srr.vector)){
    counts.subset <- read.delim(paste0(SALMON_OUTPUT.DIR, "/",srr.vector[j],"/quant.sf"))
    counts[,i] <- counts[,i] + counts.subset[match(rownames(counts),counts.subset[,1]),5]
  }
}

colSums(counts)
```

```{R, eval = F}
   PECHA_ISE6_mRNA1    PECHA_ISE6_mRNA2    PECHA_ISE6_mRNA3
           363096.0            594353.0            533546.0
PECHA_ISE6_HF_mRNA1 PECHA_ISE6_HF_mRNA2 PECHA_ISE6_HF_mRNA3
           592531.2            426703.1           1029629.3
```

### Create TPM data frame

```{R}
genelength <- read.delim(paste0(SALMON_OUTPUT.DIR, "/", groups[1,1],"/quant.sf"))
genelength <- genelength[match(rownames(counts),genelength[,1]),2]

tpm <- counts
for(i in 1:ncol(tpm)){
  tpm[,i] <- tpm[,i]/genelength
  tpm[,i] <- tpm[,i]/(sum(tpm[,i])/1000000)
}

dim(tpm)
```

```{R, eval = F}
[1] 20486     6
```

### Set group levels
```{R}
groups <- unique(groups[2:5])
groups[,1] <- factor(groups[,1], levels = groups[,1])
groups[,2] <- factor(groups[,2], levels = unique(groups[,2]))
groups[,3] <- factor(groups[,3], levels = unique(groups[,3]))
groups[,4] <- factor(groups[,4], levels = unique(groups[,4]))
```

### Identify differentially expressed genes longitudinally

edgeR and DESeq2 are both run with a FDR cutoff of <0.05 and a minimum CPM cutoff of 5 reads in the lowest sequenced sample in the data set.  

```{R}
FDRcutoff <- 0.05
cpm.cutoff <- 5/min(colSums(counts)) * 1000000

y <- DGEList(counts = counts, group = groups[,4])
y <- calcNormFactors(y)
keep <- rowSums(cpm(y) >= cpm.cutoff) >= min(table(groups[,4]))
keep.df <- as.data.frame(table(keep))
print(paste0(keep.df[keep.df[,1] == F,2]," genes excluded with CPM cutoff"))

y <- y[keep, , keep.lib.sizes = F]
design <- model.matrix(~groups[,4])
y <- estimateDisp(y , design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2:ncol(fit))

FDR <- as.data.frame(p.adjust(qlf$table$PValue, method="BH"))
rownames(FDR) <- rownames(qlf$table)
FDR.genes <- rownames(FDR[FDR[,1] < FDRcutoff, , drop = F])
print(paste0(length(FDR.genes)," DE genes identified using edgeR pairwise" ))

print(paste0(length(which(qlf$table[rownames(qlf$table) %in% FDR.genes,1] > 0) == T)," genes significantly upregulated in ",groups[1,4]))
print(paste0(length(which(qlf$table[rownames(qlf$table) %in% FDR.genes,1] < 0) == T)," genes significantly downregulated in ",groups[1,4]))

qlf$table <- as.data.frame(cbind(qlf$table,
                                 FDR))
colnames(qlf$table)[5] <- "FDR"
write.table(qlf$table[qlf$table$FDR < 0.05,],
            paste0(WORKING.DIR,"/tick_de_logFC.tsv"),
            col.names = T,
            row.names = T,
            quote = F,
            sep = "\t")
```

```{R, eval = F}
[1] "13191 genes excluded with CPM cutoff"
[1] "1793 DE genes identified using edgeR pairwise"
[1] "975 genes significantly upregulated in HF"
[1] "818 genes significantly downregulated in HF"
```

### Plot a heatmap of the up- and down-regulated genes

#### Create z-score log2TPM legend

```{R, fig,height = 2, fig.width = 7}
hmcol <- colorRampPalette(c("navyblue","white","firebrick3"))(12)
hmcolor.plot <- ggplot() + 
  geom_raster(aes(x=seq(-3,3,0.5), y=seq(-3,3,0.5), fill = seq(-3,3,0.5)))+
  scale_fill_gradientn(name = "z-score log2TPM",
                       colours=hmcol,
                       breaks=c(-3,0,3))+
  theme(legend.position="bottom")+
  guides(fill = guide_colorbar(title.position = "top"))
  

heat.legend <- g_legend(hmcolor.plot)
pdf(paste0(WORKING.DIR,"/plots/tick_hfpairwise_hm_zscorelog2tpmkey.pdf"),
    height=2,
    width=7)
grid.arrange(heat.legend)
dev.off()

png(paste0(WORKING.DIR,"/plots/tick_hfpairwise_hm_zscorelog2tpmkey.png"),
    height=2,
    width=7,
    units = "in",res=300)
grid.arrange(heat.legend)
dev.off()

grid.arrange(heat.legend)
```

![image](/images/tick_hfpairwise_hm_zscorelog2tpmkey.png)

#### Plot heatmap

```{R, fig.height=7, fig.width=5}
tpm.de <- tpm[rownames(tpm) %in% FDR.genes,]

heatmap.df <- as.data.frame(t(scale(t(log2(tpm.de + 1)))))
heatmap.df[is.na(heatmap.df)] <- 0

colcol <- as.character(groups[,2])
coldendro <- pvclust(heatmap.df, method.dist="cor", method.hclust="average", nboot=100, quiet=T)

rowcol <- ifelse(qlf$table[rownames(qlf$table) %in% FDR.genes,1] > 0,"red","blue")

pdf(paste0(WORKING.DIR,"/plots/tick_hfpairwise_zscorelog2tpm_heatmap.pdf"),
    width = 5, 
    height =7)
heatmap.2(as.matrix(heatmap.df),
          col=hmcol,
          trace="none",
          labRow=vector(mode = "character", length = nrow(heatmap.df)),
          labCol=groups[,1],
          #Rowv = as.dendrogram(rowdendro),
          Colv = as.dendrogram(coldendro),
          RowSideColors=rowcol,
          ColSideColors=colcol,
          lwid=c(1,4),
          lhei = c(1,5),
          breaks = seq(-3,3,by=.5),
          symkey=F,
          dendrogram = "column")
dev.off()

png(paste0(WORKING.DIR,"/plots/tick_hfpairwise_zscorelog2tpm_heatmap.png"),
    width = 5, 
    height = 7,
    units = "in",res=300)
heatmap.2(as.matrix(heatmap.df),
          col=hmcol,
          trace="none",
          labRow=vector(mode = "character", length = nrow(heatmap.df)),
          labCol=groups[,1],
          #Rowv = as.dendrogram(rowdendro),
          Colv = as.dendrogram(coldendro),
          RowSideColors=rowcol,
          ColSideColors=colcol,
          lwid=c(1,4),
          lhei = c(1,5),
          breaks = seq(-3,3,by=.5),
          symkey=F,
          dendrogram = "column")
dev.off()

heatmap.2(as.matrix(heatmap.df),
          col=hmcol,
          trace="none",
          labRow=vector(mode = "character", length = nrow(heatmap.df)),
          labCol=groups[,1],
          #Rowv = as.dendrogram(rowdendro),
          Colv = as.dendrogram(coldendro),
          RowSideColors=rowcol,
          ColSideColors=colcol,
          lwid=c(1,4),
          lhei = c(1,5),
          breaks = seq(-3,3,by=.5),
          symkey=F,
          dendrogram = "column")
```

![image](/images/tick_hfpairwise_zscorelog2tpm_heatmap.png)

### Conduct functional term enrichment analysis on up- and down-regulated gene subsets

```{R}
geneinfo <- read.delim(GENEINFO.PATH)

terms.pairwise.up <- as.data.frame(cbind(functionaltermenrichment(rownames(qlf$table)[qlf$table[,1] > 0],geneinfo),
                                         "upregulated"))
terms.pairwise.down <- as.data.frame(cbind(functionaltermenrichment(rownames(qlf$table)[qlf$table[,1] < 0],geneinfo),
                                     "downregulated"))

terms.colnames <- c("term","clusteroccurences","genomeoccurences","pvalue","correctedpvalue","oddsratio","group")
colnames(terms.pairwise.up) <- terms.colnames
colnames(terms.pairwise.down) <- terms.colnames
terms.pairwise <- as.data.frame(rbind(terms.pairwise.up,
                                   terms.pairwise.down))

for(i in 2:6){
  terms.pairwise[,i] <- as.numeric(as.character(terms.pairwise[,i]))
}

terms.pairwise.sig <- terms.pairwise[terms.pairwise$correctedpvalue < 0.05 & terms.pairwise$oddsratio > 1,]

write.table(terms.pairwise,
            paste0(WORKING.DIR,"/tick_hfpairwise_functionaltermenrichment.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")

write.table(terms.pairwise.sig,
            paste0(WORKING.DIR,"/tick_hfpairwise_functionaltermenrichment_sig.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")
```

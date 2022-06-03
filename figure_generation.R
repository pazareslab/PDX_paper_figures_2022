library(tidyverse)
library(progeny) 
library(GenVisR)
library(DESeq2)
library(pheatmap)
library(reshape2)

# Set the path where de .rds files are located:
path_to_rds <- "~/Desktop/PDX_paper_scripts/"

# Load clinical data:
Clinical <- readRDS(paste0(path_to_rds,"Clinical_data.rds"))
ADCnames <- Clinical %>% filter(Histology=="ADC") %>% pull(Sample)
SCCnames <- Clinical %>% filter(Histology=="SCC") %>% pull(Sample)

# Load CNV data and create PDF plots:
CNVs <- readRDS(paste0(path_to_rds,"CNVs.rds"))

ADC_NOX <- CNVs %>%
  filter(Sample %in% ADCnames & weight>1000 & cn!=2) %>% 
  select(sample=Sample, chromosome, start, end, segmean=cn)
SCC_NOX <- CNVs %>%
  filter(Sample %in% SCCnames & weight>1000 & cn!=2) %>% 
  select(sample=Sample, chromosome, start, end, segmean=cn)

pdf(paste0(path_to_rds,"ADC_CNV_PDX.pdf"),width = 20,height = 10)
cnSpec(ADC_NOX, genome="hg38",CNscale = "absolute")
dev.off()
pdf(paste0(path_to_rds,"SCC_CNV_PDX.pdf"),width = 20,height = 10)
cnSpec(SCC_NOX, genome="hg38",CNscale = "absolute")
dev.off()

# Load small variants:
Muts <- readRDS(paste0(path_to_rds,"Small_variants.rds"))

# Get info about the TMB and the top mutated genes pero histology, useful for the waterfall plot and the expression heatmap annotations:
TMB_PDL1_ADC <- Clinical %>% filter(Histology=="ADC") %>% select(Sample,"TMB High" = `TMB Consensus Status`,"PDL1 Positive" ="PDL1 Status")  %>% column_to_rownames("Sample")
top_smallvars_ADC <- Muts %>% filter(Histology == "ADC" &`Variant Category`=="Small Variant" ) %>% select(Sample,Gene) %>% group_by(Gene) %>% 
  mutate(occ=n()) %>% ungroup %>% filter(occ>=3) %>% unique %>% mutate(matrix=1) %>% select(-occ) %>% spread(Sample,matrix,fill=0) %>% column_to_rownames("Gene") %>% t %>% as.data.frame()
colnames(top_smallvars_ADC) <- top_smallvars_ADC %>% colnames %>% paste0(.," Mut")

top_CNVs_ADC <- Muts %>% filter(Histology == "ADC" &`Variant Category`=="CNV" ) %>% select(Sample,Gene) %>% group_by(Gene) %>% 
  mutate(occ=n()) %>% ungroup %>% filter(occ>=4) %>% unique %>% mutate(matrix=1) %>% select(-occ) %>% spread(Sample,matrix,fill=0) %>% column_to_rownames("Gene") %>% t %>% as.data.frame 
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}
add_to_cnvs <- outersect(ADCnames,top_CNVs_ADC %>% rownames)
for (i in 1:length(add_to_cnvs)){
  top_CNVs_ADC[add_to_cnvs[i],] <- 0
}
colnames(top_CNVs_ADC) <- top_CNVs_ADC %>% colnames %>% paste0(.," Amp")

annotation_progeny_ADC <- cbind(top_smallvars_ADC,top_CNVs_ADC,TMB_PDL1_ADC)

TMB_PDL1_SCC <- Clinical %>% filter(Histology=="SCC") %>% select(Sample,"TMB High" = `TMB Consensus Status`,"PDL1 Positive" ="PDL1 Status") %>% column_to_rownames("Sample")

top_smallvars_SCC <- Muts %>% filter(Histology == "SCC" &`Variant Category`=="Small Variant" ) %>% select(Sample,Gene) %>% group_by(Gene) %>% 
  mutate(occ=n()) %>% ungroup %>% filter(occ>=4) %>% unique %>% mutate(matrix=1) %>% select(-occ) %>% spread(Sample,matrix,fill=0) %>% column_to_rownames("Gene") %>% t %>% as.data.frame()

colnames(top_smallvars_SCC) <- top_smallvars_SCC %>% colnames %>% paste0(.," Mut")
top_CNVs_SCC <- Muts %>% filter(Histology == "SCC" &`Variant Category`=="CNV" ) %>% select(Sample,Gene) %>% group_by(Gene) %>% 
  mutate(occ=n()) %>% ungroup %>% filter(occ>=5) %>% unique %>% mutate(matrix=1) %>% select(-occ) %>% spread(Sample,matrix,fill=0) %>% column_to_rownames("Gene") %>% t %>% as.data.frame 
add_to_cnvs <- outersect(SCCnames,top_CNVs_SCC %>% rownames)
add_to_cnvs  
for (i in 1:length(add_to_cnvs)){
  top_CNVs_SCC[add_to_cnvs[i],] <- 0
}
colnames(top_CNVs_SCC) <- top_CNVs_SCC %>% colnames %>% paste0(.," Amp")
annotation_progeny_SCC <- cbind(top_smallvars_SCC,top_CNVs_SCC,TMB_PDL1_SCC)

# Load RNA-seq gene quantification data:
Exp <- readRDS(paste0(path_to_rds,"PDX_expression_readcounts.rds"))

# Per-histology filtering of low expressed and low variable genes:
vst_ADC <- Exp  %>% select(Gene,ADCnames) %>% column_to_rownames("Gene") %>%
  filter_at(.,vars(contains("TP")),any_vars(. !=0)) %>%
  filter_at(.,vars(contains("TP")),any_vars(. > 10)) %>%
  mutate(media=dplyr::select(.,contains("TP")) %>% rowMeans) %>%
  filter(media>10) %>% dplyr::select(-media) %>%  as.matrix %>%
  varianceStabilizingTransformation() %>% as.data.frame %>% mutate(variance=apply(.,1,var)) 
vst_ADC <-  vst_ADC[vst_ADC$variance >= quantile(vst_ADC$variance, c(.50)), ] %>%
  select(-variance) 

vst_SCC <- Exp  %>% select(Gene,SCCnames) %>% column_to_rownames("Gene") %>%
  filter_at(.,vars(contains("TP")),any_vars(. !=0)) %>%
  filter_at(.,vars(contains("TP")),any_vars(. > 10)) %>%
  mutate(media=dplyr::select(.,contains("TP")) %>% rowMeans) %>%
  filter(media>10) %>% dplyr::select(-media) %>%  as.matrix %>%
  varianceStabilizingTransformation() %>% as.data.frame %>% mutate(variance=apply(.,1,var)) 
vst_SCC <-  vst_SCC[vst_SCC$variance >= quantile(vst_SCC$variance, c(.50)), ] %>%
  select(-variance) 

# Run Progeny (pathway activity analysis):
prog_ADC <- vst_ADC %>% as.matrix %>%   progeny() %>% t 
prog_SCC <- vst_SCC %>% as.matrix %>%  progeny() %>% t 

# Prepare heatmap annotations:
common_annotations <- intersect(annotation_progeny_ADC %>% colnames(), annotation_progeny_SCC %>% colnames())
unique_annotations <- outersect(annotation_progeny_ADC %>% colnames(), annotation_progeny_SCC %>% colnames())

annotation_colors_ADC <-  list("TP53 Mut" = c("1"="black", "0" = "white"),
                               "MYC Amp" = c("1"="darkred", "0" = "white"),
                               "TMB High" = c("1"="darkblue", "0" = "white"),
                               "PDL1 Positive" = c("1" = "darkgreen", "0" = "white"),
                               "EGFR Mut" = c("1" = "#B15928", "0" = "white"),                             
                               "KRAS Mut" = c("1" = "#FFFF33", "0" = "white"),                             
                               "MGA Mut" = c("1" = "#FF7F00", "0" = "white"),                             
                               "PIK3CA Mut" = c("1" = "#984EA3", "0" = "white"),                             
                               "RB1 Mut" = c("1" = "#DECBE4", "0" = "white"),                             
                               "RET Mut" = c("1" = "#FFD92F", "0" = "white"),                             
                               "STK11 Mut" = c("1" = "#7FC97F", "0" = "white"),                             
                               "FGFR1 Amp" = c("1" = "#F781BF", "0" = "white"))                             

rownames_adc <- annotation_progeny_ADC %>% rownames
annotation_progeny_ADC <- annotation_progeny_ADC %>% apply(.,2,as.character) %>% as.data.frame
rownames(annotation_progeny_ADC) <- rownames_adc
pdf(paste0(path_to_rds,"ADC_progeny_PDX.pdf"), width = 10, height = 8)
prog_ADC %>% pheatmap(annotation_colors = annotation_colors_ADC,annotation_col = annotation_progeny_ADC, clustering_distance_cols = "correlation",annotation_legend = F) # ESTE!!! (4 grupos: 1a)PIK3CA+KMT2D 1b)PIK3CA 2a)Nada especial? 2b) MYC+CDKN2A
dev.off()

annotation_colors_SCC <-  list("TP53 Mut" = c("1"="black", "0" = "white"),
                               "MYC Amp" = c("1"="darkred", "0" = "white"),
                               "TMB High" = c("1"="darkblue", "0" = "white"),
                               "PDL1 Positive" = c("1" = "darkgreen", "0" = "white"),
                               "PIK3CA Amp" = c("1" = "#7FC97F", "0" = "white"),
                               "MET Amp" = c("1" = "#DECBE4", "0" = "white"),
                               "NOTCH1 Mut" = c("1" = "#E5D8BD", "0" = "white"),
                               "KMT2D Mut" = c("1" = "#A65628", "0" = "white"),
                               "KEAP1 Mut" = c("1" = "#FFFF33", "0" = "white"),
                               "GNAS Mut" = c("1" = "#B15928", "0" = "white"),
                               "FAT1 Mut" = c("1" = "#984EA3", "0" = "white"),
                               "ESR1 Mut" = c("1" = "#FFD92F", "0" = "white"),
                               "CDKN2A Mut" = c("1" = "#FF7F00", "0" = "white"),
                               "ARID1A Mut" = c("1" = "#F781BF", "0" = "white"))
rownames_scc <- annotation_progeny_SCC %>% rownames
annotation_progeny_SCC <- annotation_progeny_SCC %>% apply(.,2,as.character) %>% as.data.frame
rownames(annotation_progeny_SCC) <- rownames_scc
pdf(paste0(path_to_rds,"SCC_progeny_PDX.pdf"), width = 10, height = 8)
prog_SCC %>% pheatmap(annotation_colors = annotation_colors_SCC,annotation_col = annotation_progeny_SCC, clustering_distance_cols = "correlation",annotation_legend = F) # ESTE!!! (4 grupos: 1a)PIK3CA+KMT2D 1b)PIK3CA 2a)Nada especial? 2b) MYC+CDKN2A
dev.off()

# Get Waterfall mutation plots per histology:

genelist <- readRDS(paste0(path_to_rds,"gene_list.rds"))
oncoADC <-   Muts %>% filter(Histology=="ADC") %>% 
  select(sample=Sample, gene=Gene, variant_type=`Extra info`, variant_class=`Variant effect`) %>% as.data.frame 
most_deleterious <- c("Double","Fusion","Startloss","Stopgain","Frameshift","Splicing","Loss","Missense","Gain","Nonframeshift")
colors <- c("black","slategrey","#762A83","magenta","darkorange","dodgerblue","blue","darkgreen","darkred","yellowgreen")
clinicalADC <- Clinical %>% filter(Histology=="ADC") %>% select(sample=Sample, "TMB" = "TMB Consensus Status") #%>% gather(key="variable", value= value )
clinicalADC <- melt(clinicalADC, id.vars = c("sample"))
clinicalADC <- clinicalADC %>% mutate(value = replace(value, value==0, "Low (<10)")) %>% 
  mutate(value = replace(value, value==1, "High (>=10)"))

pdf(paste0(path_to_rds,"waterfall_ADC_PDX.pdf"), width = 30, height = 15)
waterfall(oncoADC,
          plotGenes = genelist,
          mainDropMut = TRUE,
          rmvSilent = TRUE,
          mainLabelCol = "variant_type",
          clinData = clinicalADC,
          fileType = "Custom",
          variant_class_order = most_deleterious,
          clinVarCol = c("High (>=10)" =  "darkblue",
                         "Low (<10)" = "lightblue"),          
          mainPalette = colors,
          mainXlabel = TRUE,
          section_heights = c(0,1,0.05),
          main_geneLabSize = 15,
          plotMutBurden = FALSE,
          mainLabelSize = 5,
          mainGrid = F
)
dev.off()

oncoSCC <-   Muts %>% filter(Sample %in% SCCnames) %>% 
  select(sample=Sample, gene=Gene, variant_type=`Extra info`, variant_class=`Variant effect`) %>% as.data.frame
clinicalSCC <- Clinical %>%filter(Sample %in% SCCnames)%>% select(sample=Sample, "TMB" = "TMB Consensus Status", Histology) %>% gather(key="variable", value= value,-sample )
clinicalSCC<- clinicalSCC %>% mutate(value = replace(value, value==0, "Low (<10)")) %>% 
  mutate(value = replace(value, value==1, "High (>=10)"))
pdf(paste0(path_to_rds,"waterfall_SCC_PDX.pdf"), width = 30, height = 15)
waterfall(oncoSCC,
          plotGenes = genelist,
          mainDropMut = TRUE,
          rmvSilent = TRUE,
          mainLabelCol = "variant_type",
          clinData = clinicalSCC,
          fileType = "Custom",
          variant_class_order = most_deleterious,
          clinVarCol = c("High (>=10)" =  "darkblue",
                         "Low (<10)" = "lightblue"),          
          mainPalette = colors,
          mainXlabel = TRUE,
          section_heights = c(0,1,0.05),
          main_geneLabSize = 15,
          plotMutBurden = FALSE,
          mainLabelSize = 4.5, #4,5
          mainGrid = F
)

dev.off()



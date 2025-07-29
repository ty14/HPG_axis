library(RRHO2)
library(tidyverse)
#Sources 
source("females/analysis/06_mostvar_genes_allTissues.R")
ARC <- arc_var %>% rownames_to_column(., var = "genes")
VMN <- vmn_var %>% rownames_to_column(., var = "genes")
PIT <- pit_var %>% rownames_to_column(., var = "genes")
OVT <- ovt_var %>% rownames_to_column(., var = "genes")

OVTv <- OVT[OVT$genes %in% VMN$genes,] #2406
OVTa <- OVT[OVT$genes %in% ARC$genes,] #2527
OVTp <- OVT[OVT$genes %in% PIT$genes,] #2755

VMN <- VMN[VMN$genes %in% OVTv$genes,]
ARC <- ARC[ARC$genes %in% OVTa$genes,]
PIT <- PIT[PIT$genes %in% OVTp$genes,]


#first vmn
fe_vmn <- readRDS("females/clean_data/limma_eFDR_VMN.RDS") %>% filter(genes %in% VMN$genes)
fe_ovt <- readRDS("females/clean_data/limma_eFDR_OVT.RDS") %>% filter(genes %in% OVTv$genes)
ovtx <- fe_ovt %>% mutate(gene = row_number())
vmnx <- fe_vmn %>% mutate(gene = row_number())



vmnx$gene <- paste0("Gene", vmnx$gene)
head(vmnx)
ovtx$gene <- paste0("Gene", ovtx$gene)
head(ovtx)
#List 1
Gene = ovtx$gene

list1_pvalue_1_200 <- ovtx %>% filter(logFC > 0.2 & P.Value < 0.05) 
list1_pvalue_1_200 <-list1_pvalue_1_200$P.Value## up-regulated genes
list1_pvalue_201_400 <- ovtx %>% filter(logFC < 0.2 & P.Value < 0.05)
list1_pvalue_201_400 <-list1_pvalue_201_400$P.Value ## down-regulated genes
an <- subset(ovtx, !(P.Value %in% list1_pvalue_1_200))
an <- subset(an, !(P.Value %in% list1_pvalue_201_400))
list1_pvalue_401_2000 <- an$P.Value  ## non-changed genes
list1_DDE <- c(-log10(list1_pvalue_1_200), -log10(list1_pvalue_201_400) * (-1), -log10(list1_pvalue_401_2000) * sample(c(1,-1), length(list1_pvalue_401_2000), replace = TRUE))


lapply(list1_DDE, head)
gene_list1 <- data.frame(Genes=Gene,DDE = list1_DDE, stringsAsFactors = FALSE)

#List2
list2_pvalue_1_200 <-vmnx %>% filter(logFC > 0.2 & P.Value < 0.05) 
list2_pvalue_1_200 <-list2_pvalue_1_200$P.Value## up-regulated genes
list2_pvalue_201_400 <- vmnx %>% filter(logFC < 0.2 & P.Value < 0.05)
list2_pvalue_201_400 <-list2_pvalue_201_400$P.Value ## down-regulated genes
pn <- subset(vmnx, !(P.Value %in% list2_pvalue_1_200))
pn <- subset(pn, !(P.Value %in% list2_pvalue_201_400))
list2_pvalue_401_2000 <- pn$P.Value  ## non-changed genes
list2_DDE <- c(-log10(list2_pvalue_1_200), -log10(list2_pvalue_201_400) * (-1), -log10(list2_pvalue_401_2000) * sample(c(1,-1), length(list2_pvalue_401_2000), replace = TRUE))



gene_list2 <- data.frame(Genes=Gene,DDE = list2_DDE, stringsAsFactors = FALSE)

library(RRHO2)
RRHO_obj <-  RRHO2_initialize(gene_list1, gene_list2, labels = c("OVT", "VMN"), log10.ind=TRUE)


RRHO2_heatmap(RRHO_obj)

# -----------------------------------

#first arc
fe_arc <- readRDS("females/clean_data/limma_eFDR_ARC.RDS") %>% filter(genes %in% ARC$genes)
fe_ovt <- readRDS("females/clean_data/limma_eFDR_OVT.RDS") %>% filter(genes %in% OVTa$genes)
arcx <- fe_arc %>% mutate(gene = row_number())
ovtx <- fe_ovt %>% mutate(gene = row_number())



arcx$gene <- paste0("Gene", arcx$gene)
head(arcx)
ovtx$gene <- paste0("Gene", ovtx$gene)
head(ovtx)
#List 1
Gene = ovtx$gene

list1_pvalue_1_200 <- ovtx %>% filter(logFC > 0.2 & P.Value < 0.05) 
list1_pvalue_1_200 <-list1_pvalue_1_200$P.Value## up-regulated genes
list1_pvalue_201_400 <- ovtx %>% filter(logFC < 0.2 & P.Value < 0.05)
list1_pvalue_201_400 <-list1_pvalue_201_400$P.Value ## down-regulated genes
an <- subset(ovtx, !(P.Value %in% list1_pvalue_1_200))
an <- subset(an, !(P.Value %in% list1_pvalue_201_400))
list1_pvalue_401_2000 <- an$P.Value  ## non-changed genes
list1_DDE <- c(-log10(list1_pvalue_1_200), -log10(list1_pvalue_201_400) * (-1), -log10(list1_pvalue_401_2000) * sample(c(1,-1), length(list1_pvalue_401_2000), replace = TRUE))


lapply(list1_DDE, head)
gene_list1 <- data.frame(Genes=Gene,DDE = list1_DDE, stringsAsFactors = FALSE)

#List2
list2_pvalue_1_200 <-arcx %>% filter(logFC > 0.2 & P.Value < 0.05) 
list2_pvalue_1_200 <-list2_pvalue_1_200$P.Value## up-regulated genes
list2_pvalue_201_400 <- arcx %>% filter(logFC < 0.2 & P.Value < 0.05)
list2_pvalue_201_400 <-list2_pvalue_201_400$P.Value ## down-regulated genes
pn <- subset(arcx, !(P.Value %in% list2_pvalue_1_200))
pn <- subset(pn, !(P.Value %in% list2_pvalue_201_400))
list2_pvalue_401_2000 <- pn$P.Value  ## non-changed genes
list2_DDE <- c(-log10(list2_pvalue_1_200), -log10(list2_pvalue_201_400) * (-1), -log10(list2_pvalue_401_2000) * sample(c(1,-1), length(list2_pvalue_401_2000), replace = TRUE))



gene_list2 <- data.frame(Genes=Gene,DDE = list2_DDE, stringsAsFactors = FALSE)

library(RRHO2)
RRHO_obj <-  RRHO2_initialize(gene_list1, gene_list2, labels = c("OVT", "ARC"), log10.ind=TRUE)


RRHO2_heatmap(RRHO_obj)


#first pit
fe_pit <- readRDS("females/clean_data/limma_eFDR_PIT.RDS") %>% filter(genes %in% PIT$genes)
fe_ovt <- readRDS("females/clean_data/limma_eFDR_OVT.RDS") %>% filter(genes %in% OVTp$genes)
pitx <- fe_pit %>% mutate(gene = row_number())
ovtx <- fe_ovt %>% mutate(gene = row_number())



pitx$gene <- paste0("Gene", pitx$gene)
head(pitx)
ovtx$gene <- paste0("Gene", ovtx$gene)
head(ovtx)
#List 1
Gene = ovtx$gene

list1_pvalue_1_200 <- ovtx %>% filter(logFC > 0.2 & P.Value < 0.05) 
list1_pvalue_1_200 <-list1_pvalue_1_200$P.Value## up-regulated genes
list1_pvalue_201_400 <- ovtx %>% filter(logFC < 0.2 & P.Value < 0.05)
list1_pvalue_201_400 <-list1_pvalue_201_400$P.Value ## down-regulated genes
an <- subset(ovtx, !(P.Value %in% list1_pvalue_1_200))
an <- subset(an, !(P.Value %in% list1_pvalue_201_400))
list1_pvalue_401_2000 <- an$P.Value  ## non-changed genes
list1_DDE <- c(-log10(list1_pvalue_1_200), -log10(list1_pvalue_201_400) * (-1), -log10(list1_pvalue_401_2000) * sample(c(1,-1), length(list1_pvalue_401_2000), replace = TRUE))


lapply(list1_DDE, head)
gene_list1 <- data.frame(Genes=Gene,DDE = list1_DDE, stringsAsFactors = FALSE)

#List2
list2_pvalue_1_200 <-pitx %>% filter(logFC > 0.2 & P.Value < 0.05) 
list2_pvalue_1_200 <-list2_pvalue_1_200$P.Value## up-regulated genes
list2_pvalue_201_400 <- pitx %>% filter(logFC < 0.2 & P.Value < 0.05)
list2_pvalue_201_400 <-list2_pvalue_201_400$P.Value ## down-regulated genes
pn <- subset(pitx, !(P.Value %in% list2_pvalue_1_200))
pn <- subset(pn, !(P.Value %in% list2_pvalue_201_400))
list2_pvalue_401_2000 <- pn$P.Value  ## non-changed genes
list2_DDE <- c(-log10(list2_pvalue_1_200), -log10(list2_pvalue_201_400) * (-1), -log10(list2_pvalue_401_2000) * sample(c(1,-1), length(list2_pvalue_401_2000), replace = TRUE))



gene_list2 <- data.frame(Genes=Gene,DDE = list2_DDE, stringsAsFactors = FALSE)

library(RRHO2)
RRHO_obj <-  RRHO2_initialize(gene_list1, gene_list2, labels = c("OVT", "PIT"), log10.ind=TRUE)


RRHO2_heatmap(RRHO_obj)




#first vmn and arc 

VMNa <- VMN[VMN$genes %in% ARC$genes,]
ARCv<- ARC[ARC$genes %in% VMN$genes,]
fe_vmn <- readRDS("females/clean_data/limma_eFDR_VMN.RDS") %>% filter(genes %in% VMNa$genes)
fe_ovt <- readRDS("females/clean_data/limma_eFDR_ARC.RDS") %>% filter(genes %in% ARCv$genes)

ovtx <- fe_ovt %>% mutate(gene = row_number())
vmnx <- fe_vmn %>% mutate(gene = row_number())



vmnx$gene <- paste0("Gene", vmnx$gene)
head(vmnx)
ovtx$gene <- paste0("Gene", ovtx$gene)
head(ovtx)
#List 1
Gene = ovtx$gene

list1_pvalue_1_200 <- ovtx %>% filter(logFC > 0.2 & P.Value < 0.05) 
list1_pvalue_1_200 <-list1_pvalue_1_200$P.Value## up-regulated genes
list1_pvalue_201_400 <- ovtx %>% filter(logFC < 0.2 & P.Value < 0.05)
list1_pvalue_201_400 <-list1_pvalue_201_400$P.Value ## down-regulated genes
an <- subset(ovtx, !(P.Value %in% list1_pvalue_1_200))
an <- subset(an, !(P.Value %in% list1_pvalue_201_400))
list1_pvalue_401_2000 <- an$P.Value  ## non-changed genes
list1_DDE <- c(-log10(list1_pvalue_1_200), -log10(list1_pvalue_201_400) * (-1), -log10(list1_pvalue_401_2000) * sample(c(1,-1), length(list1_pvalue_401_2000), replace = TRUE))


lapply(list1_DDE, head)
gene_list1 <- data.frame(Genes=Gene,DDE = list1_DDE, stringsAsFactors = FALSE)

#List2
list2_pvalue_1_200 <-vmnx %>% filter(logFC > 0.2 & P.Value < 0.05) 
list2_pvalue_1_200 <-list2_pvalue_1_200$P.Value## up-regulated genes
list2_pvalue_201_400 <- vmnx %>% filter(logFC < 0.2 & P.Value < 0.05)
list2_pvalue_201_400 <-list2_pvalue_201_400$P.Value ## down-regulated genes
pn <- subset(vmnx, !(P.Value %in% list2_pvalue_1_200))
pn <- subset(pn, !(P.Value %in% list2_pvalue_201_400))
list2_pvalue_401_2000 <- pn$P.Value  ## non-changed genes
list2_DDE <- c(-log10(list2_pvalue_1_200), -log10(list2_pvalue_201_400) * (-1), -log10(list2_pvalue_401_2000) * sample(c(1,-1), length(list2_pvalue_401_2000), replace = TRUE))



gene_list2 <- data.frame(Genes=Gene,DDE = list2_DDE, stringsAsFactors = FALSE)

library(RRHO2)
RRHO_obj <-  RRHO2_initialize(gene_list1, gene_list2, labels = c("ARC", "VMN"), log10.ind=TRUE)


RRHO2_heatmap(RRHO_obj)

# -----------------------------------

#first arc
PITa <- PIT[PIT$genes %in% ARC$genes,]
ARCp<- ARC[ARC$genes %in% PIT$genes,]

fe_pit <- readRDS("females/clean_data/limma_eFDR_PIT.RDS") %>% filter(genes %in% PITa$genes)
fe_ovt <- readRDS("females/clean_data/limma_eFDR_ARC.RDS") %>% filter(genes %in% ARCp$genes)
arcx <- fe_pit %>% mutate(gene = row_number())
ovtx <- fe_ovt %>% mutate(gene = row_number())



arcx$gene <- paste0("Gene", arcx$gene)
head(arcx)
ovtx$gene <- paste0("Gene", ovtx$gene)
head(ovtx)
#List 1
Gene = ovtx$gene

list1_pvalue_1_200 <- ovtx %>% filter(logFC > 0.2 & P.Value < 0.05) 
list1_pvalue_1_200 <-list1_pvalue_1_200$P.Value## up-regulated genes
list1_pvalue_201_400 <- ovtx %>% filter(logFC < 0.2 & P.Value < 0.05)
list1_pvalue_201_400 <-list1_pvalue_201_400$P.Value ## down-regulated genes
an <- subset(ovtx, !(P.Value %in% list1_pvalue_1_200))
an <- subset(an, !(P.Value %in% list1_pvalue_201_400))
list1_pvalue_401_2000 <- an$P.Value  ## non-changed genes
list1_DDE <- c(-log10(list1_pvalue_1_200), -log10(list1_pvalue_201_400) * (-1), -log10(list1_pvalue_401_2000) * sample(c(1,-1), length(list1_pvalue_401_2000), replace = TRUE))


lapply(list1_DDE, head)
gene_list1 <- data.frame(Genes=Gene,DDE = list1_DDE, stringsAsFactors = FALSE)

#List2
list2_pvalue_1_200 <-arcx %>% filter(logFC > 0.2 & P.Value < 0.05) 
list2_pvalue_1_200 <-list2_pvalue_1_200$P.Value## up-regulated genes
list2_pvalue_201_400 <- arcx %>% filter(logFC < 0.2 & P.Value < 0.05)
list2_pvalue_201_400 <-list2_pvalue_201_400$P.Value ## down-regulated genes
pn <- subset(arcx, !(P.Value %in% list2_pvalue_1_200))
pn <- subset(pn, !(P.Value %in% list2_pvalue_201_400))
list2_pvalue_401_2000 <- pn$P.Value  ## non-changed genes
list2_DDE <- c(-log10(list2_pvalue_1_200), -log10(list2_pvalue_201_400) * (-1), -log10(list2_pvalue_401_2000) * sample(c(1,-1), length(list2_pvalue_401_2000), replace = TRUE))



gene_list2 <- data.frame(Genes=Gene,DDE = list2_DDE, stringsAsFactors = FALSE)

library(RRHO2)
RRHO_obj <-  RRHO2_initialize(gene_list1, gene_list2, labels = c("ARC", "PIT"), log10.ind=TRUE)


RRHO2_heatmap(RRHO_obj)




#first arc
PITv <- PIT[PIT$genes %in% VMN$genes,]
VMNp<- VMN[VMN$genes %in% PIT$genes,]

fe_pit <- readRDS("females/clean_data/limma_eFDR_PIT.RDS") %>% filter(genes %in% PITv$genes)
fe_ovt <- readRDS("females/clean_data/limma_eFDR_VMN.RDS") %>% filter(genes %in% VMNp$genes)
arcx <- fe_pit %>% mutate(gene = row_number())
ovtx <- fe_ovt %>% mutate(gene = row_number())



arcx$gene <- paste0("Gene", arcx$gene)
head(arcx)
ovtx$gene <- paste0("Gene", ovtx$gene)
head(ovtx)
#List 1
Gene = ovtx$gene

list1_pvalue_1_200 <- ovtx %>% filter(logFC > 0.2 & P.Value < 0.05) 
list1_pvalue_1_200 <-list1_pvalue_1_200$P.Value## up-regulated genes
list1_pvalue_201_400 <- ovtx %>% filter(logFC < 0.2 & P.Value < 0.05)
list1_pvalue_201_400 <-list1_pvalue_201_400$P.Value ## down-regulated genes
an <- subset(ovtx, !(P.Value %in% list1_pvalue_1_200))
an <- subset(an, !(P.Value %in% list1_pvalue_201_400))
list1_pvalue_401_2000 <- an$P.Value  ## non-changed genes
list1_DDE <- c(-log10(list1_pvalue_1_200), -log10(list1_pvalue_201_400) * (-1), -log10(list1_pvalue_401_2000) * sample(c(1,-1), length(list1_pvalue_401_2000), replace = TRUE))


lapply(list1_DDE, head)
gene_list1 <- data.frame(Genes=Gene,DDE = list1_DDE, stringsAsFactors = FALSE)

#List2
list2_pvalue_1_200 <-arcx %>% filter(logFC > 0.2 & P.Value < 0.05) 
list2_pvalue_1_200 <-list2_pvalue_1_200$P.Value## up-regulated genes
list2_pvalue_201_400 <- arcx %>% filter(logFC < 0.2 & P.Value < 0.05)
list2_pvalue_201_400 <-list2_pvalue_201_400$P.Value ## down-regulated genes
pn <- subset(arcx, !(P.Value %in% list2_pvalue_1_200))
pn <- subset(pn, !(P.Value %in% list2_pvalue_201_400))
list2_pvalue_401_2000 <- pn$P.Value  ## non-changed genes
list2_DDE <- c(-log10(list2_pvalue_1_200), -log10(list2_pvalue_201_400) * (-1), -log10(list2_pvalue_401_2000) * sample(c(1,-1), length(list2_pvalue_401_2000), replace = TRUE))



gene_list2 <- data.frame(Genes=Gene,DDE = list2_DDE, stringsAsFactors = FALSE)

library(RRHO2)
RRHO_obj <-  RRHO2_initialize(gene_list1, gene_list2, labels = c("VMN", "PIT"), log10.ind=TRUE)


RRHO2_heatmap(RRHO_obj)

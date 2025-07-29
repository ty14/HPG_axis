library(limma)
library(DESeq2)
library(edgeR)
library(AnnotationDbi)
library(tidyverse)
organism = 'org.Rn.eg.db'
library('org.Rn.eg.db')
library(limma)
library(DESeq2)
library(edgeR)
library(AnnotationDbi)
library(tidyverse)
organism = 'org.Rn.eg.db'
library('org.Rn.eg.db')


# remove NMX03B 
df <- read.table("raw_data/all_gonad_counts_Star.gff", col.names = c("symbol", "count","id"))
head(df)
dfx <- df %>% pivot_wider(names_from = id, values_from = count) %>% column_to_rownames(., var = "symbol")
colnames(dfx) <- substr(colnames(dfx), 7, 12)
# dfxx <- dfx %>% rownames_to_column(., var = "symbol")
# saveRDS(dfxx, "raw_data/gonads_clean_counts.RDS")

#counts
g <- dfx

# sample information 
id <- read_csv("raw_data/traitData.csv")
colnames(id)
idx <- id %>% dplyr::select(1:4,15)

#coldata vmn
g_data <- idx %>% filter(SampleName != "NMX03B")
g_datax <- g_data %>% column_to_rownames(., var = "SampleName")
g_datax$bi_group <- ifelse(g_datax$Group == "Blue", 0, 1)

fg_datax <- g_datax %>% filter(Sex == "M") 

#checks
all(row.names(fg_datax) %in% colnames(g)) #check 
g <- g[,rownames(fg_datax)]
all(rownames(fg_datax) == colnames(g)) #check


#normalize and filter with all groups -
# I normalize with all groups since I want to compare them after

g<- g[!is.na(rowSums(g)),]
# rownames(g) <- rownames

d = apply(g, 2, as.numeric)
dim(d)
# [1] 34039 23
d0= DGEList(d, group = fg_datax$Group)
dim(d0)
rownames(d0) <- rownames(g)
d0 <- calcNormFactors(d0)

cutoff <- 10
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
#[1] 9613   23

fg_datax$Group %>%
  factor(.,levels = c("Blue", "Orange")) -> group.dl


design.dl <- model.matrix(~ 0 + group.dl)
colnames(design.dl) -> mycolnames

v.dl = voom(dge.dl, design.dl, plot = T)
vfit.dl = lmFit(v.dl, design.dl)

contr.matrix <- makeContrasts(group.dlBlue-group.dlOrange,
                              levels=design.dl)

vfit.dl2 <- contrasts.fit(vfit.dl, contr.matrix)

efit.dl2 = eBayes(vfit.dl2)
plotSA(efit.dl2, main="Final model: Mean-variance trend")

p.dl.limma2 = efit.dl2[["p.value"]]
head(p.dl.limma2)

saveRDS(v.dl, "males/clean_data/limma_vdl_TES.RDS")




# How many random sampling
R = 5000
set.seed(312)

#to store pvalues in
p.dl.rand = vector('list',length = R)

# to store "t" values (coefficients)
p.dl.rand.t = vector('list',length = R)

for( g in 1 : R){
  print(paste("Starting on Permutation", g))
  
  # Randomize the traits
  
  group.dl.rand = sample(group.dl)
  
  # Model
  design.dl.rand = model.matrix(~0 + group.dl.rand)
  colnames(design.dl.rand) <- mycolnames
  
  # Calculate p-values based on randomized traits
  v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
  vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
  
  vfit.dl.rand2 <- contrasts.fit(vfit.dl.rand, contr.matrix)
  
  efit.dl.rand2 = eBayes(vfit.dl.rand2)
  
  p.dl.rand[[g]] = efit.dl.rand2[["p.value"]]
  p.dl.rand.t[[g]] = efit.dl.rand2[["t"]]
  
}

#Counting how many observations are above/below observed values
q.dl <- Reduce(`+`, lapply(p.dl.rand, \(x) {
  (x < p.dl.limma2)
}))

q.dl

q.dl = q.dl / R
colnames(q.dl) = 'bluevs.orange'
q.dl = as.data.frame(q.dl)


efit.dl2[["p.value"]] <- q.dl
row.names(q.dl) <- NULL
sum(duplicated(row.names(efit.dl2$coefficients)))

saveRDS(q.dl,("males/clean_data/limma_vdl_eFDRcutoffs_TES.RDS"))

##### Analysis pulling genes out for each contrast 
tmp1 <- contrasts.fit(efit.dl2, coef = 1) #Blue vs. Orange


topTable(tmp1, sort.by = "P", n = Inf) %>% 
  rownames_to_column('genes') -> bvo


saveRDS(bvo,"males/clean_data/limma_eFDR_testis.RDS")

fe_bvo <- readRDS("males/clean_data/limma_eFDR_testis.RDS") 


#quick look at numbers 
fe_bvo <- readRDS("males/clean_data/limma_eFDR_testis.RDS") %>% 
  distinct(.) %>% 
  filter(.,abs(logFC) >= 0.2) %>%
  filter(.,P.Value <0.05)

#quick look at number of genes
fe_bvo %>% filter(., P.Value<0.05) %>% 
  summarise(.,Up = sum(logFC>0.2),
            Down = sum(logFC<0.2)) %>% 
  mutate(.,Total = Up + Down)

# Up Down Total
# 1    6     7

hist(fe_bvo$logFC)

##top genes 
fe_bvo_up <- fe_bvo %>% 
  filter(logFC > 0.5)

fe_bvo_up %>% arrange(-logFC) %>% head(., 20)

fe_bvo_down <- fe_bvo %>% 
  filter(logFC < 0.5)

fe_bvo_down %>% arrange(logFC) %>% head(., 20)

# 
# # go analysis
# library(clusterProfiler)
# library(AnnotationDbi)
# library(org.Rn.eg.db)
# library(tidyverse)
# 
# #Descenders 
# my_logFC_threshold = 0.5
# 
# y1a <- readRDS("clean_data/limma_eFDR_females.RDS") %>% 
#   distinct(.) %>% 
#   filter(.,abs(logFC) >= my_logFC_threshold) %>%
#   filter(.,P.Value <0.05) 
# 
# y1a
# 
# 
# y1a <- y1a %>% 
#   mutate(entrez = mapIds(org.Rn.eg.db, keys = y1a$genes,
#                          column = "ENTREZID", keytype = "SYMBOL")) %>% 
#   filter(!is.na(entrez))
# 
# source("functions/gettop10GO_Rat.R")
# 
# gettop10GO(y1a, my_showCategory) %>% 
#   mutate(tissue = "ovaries") -> top10go1
# 
# write.csv(top10go1, "clean_data/GOTerms_ovaries_limma.csv", row.names = F)



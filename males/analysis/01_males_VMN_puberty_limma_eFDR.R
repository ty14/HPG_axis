library(limma)
library(DESeq2)
library(edgeR)
library(AnnotationDbi)
library(tidyverse)
organism = 'org.Rn.eg.db'
library('org.Rn.eg.db')

# sample information 
id <- read_csv("raw_data/traitData.csv")
colnames(id)
idx <- id %>% select(1:4,15)

#counts 
v <- readRDS("brain/clean_data/vmn_counts.RDS") %>%
  select(-region) %>% 
  column_to_rownames(., "ensgene")

colnames(v) <- paste("V", colnames(v), sep = "_")

#coldata vmn
v_data <- idx
v_data$SampleName <- paste("V", v_data$SampleName, sep = "_")
v_datax <- v_data %>% column_to_rownames(., var = "SampleName")
v_datax$bi_sex <- ifelse(v_datax$Sex == "M", 0, 1)
v_datax$bi_group <- ifelse(v_datax$Group == "Blue", 0, 1)


#checks
all(row.names(v_datax) %in% colnames(v)) #check 
v_datax <- v_datax[colnames(v),]
v <- v[,rownames(v_datax)]
all(rownames(v_datax) == colnames(v)) #check



#normalize and filter with all groups -
# I normalize with all groups since I want to compare them after
d = apply(v, 2, as.numeric)
dim(d)
# [1] 34039    36
d0= DGEList(d, group = v_datax$Sex)
dim(d0)
rownames(d0) <- rownames(v)
d0 <- calcNormFactors(d0)

cutoff <- 10
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
# [1] 10755    36

#get just females: 
males <- v_datax %>% filter(Sex == "M")

males$bi_group -> group.dl

puberty.dl <- (males$Puberty-min(males$Puberty))/(max(males$Puberty)-min(males$Puberty))

dge.m <- dge.dl[, dge.dl$samples$group %in% c("M")]
dge.m$samples$group <- droplevels(dge.m$samples$group)
dge.m$samples$group
dge.dl<- dge.m
dge.dl$samples$group

all(row.names(males) %in% colnames(dge.dl$counts)) #check 


design.dl <- model.matrix(~group.dl+puberty.dl)
colnames(design.dl) -> mycolnames

v.dl = voom(dge.dl, design.dl, plot = T)
vfit.dl = lmFit(v.dl, design.dl)
efit.dl = eBayes(vfit.dl)



p.dl.limma = efit.dl[["p.value"]]


saveRDS(v.dl, "males/clean_data/limma_vdl_males_VMN_puberty.RDS")


# How many random sampling
R = 5000


p.dl.rand = vector('list',length = R)

for( g in 1 : R){
  print(paste("Starting on Permutation", g))
  
  # Randomize the traits
  group.dl.rand = sample(group.dl)
  puberty.dl.rand = sample(puberty.dl)
  # Model
  design.dl.rand = model.matrix(~group.dl.rand+puberty.dl.rand)
  colnames(design.dl.rand) <- mycolnames
  
  # Calculate p-values based on randomized traits
  v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
  vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
  
  efit.dl.rand = eBayes(vfit.dl.rand)
  
  p.dl.rand[[g]] = efit.dl.rand[["p.value"]]
  
  
}


#Counting how many observations are above/below observed values

q.dl <- Reduce(`+`, lapply(p.dl.rand, \(x) {
  (x < p.dl.limma)
}))

q.dl

q.dl = q.dl / R
colnames(q.dl) = mycolnames
q.dl = as.data.frame(q.dl)


## Replacing pvalues from limma with permutation p-values from pvalues?
efit.dl[["p.value"]] <- q.dl
row.names(q.dl) <- NULL
sum(duplicated(row.names(efit.dl$coefficients)))

saveRDS(q.dl,("males/clean_data/limma_vdl_eFDRcutoffs_males_VMN_puberty.RDS"))


tmp1 <- contrasts.fit(efit.dl, coef = 2)#group
tmp2 <- contrasts.fit(efit.dl, coef = 3)#puberty 

limma_list <- list()
topTable(tmp1, sort.by = "P", n = Inf) %>% 
  rownames_to_column('genes') -> limma_list$group

topTable(tmp2, sort.by = "P", n = Inf) %>% 
  rownames_to_column('genes') -> limma_list$puberty


saveRDS(limma_list,"males/clean_data/limma_males_VMN_puberty.RDS")


my_logFC_threshold = 0.5 ## what should this be?

limma_list<- readRDS("males/clean_data/limma_males_VMN_puberty.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) 

#colnames(limma_list$status)
limma_list %>% 
  map(~summarise(.,Up = sum(logFC>0),
                 Down = sum(logFC<0))) %>% 
  map(~mutate(.,Total = Up + Down))



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
a <- readRDS("brain/clean_data/arc_counts.RDS") %>%
  select(-region) %>% 
  column_to_rownames(., "ensgene")

colnames(a) <- paste("A", colnames(a), sep = "_")

#coldata vmn
a_data <- idx
a_data$SampleName <- paste("A", a_data$SampleName, sep = "_")
a_datax <- a_data %>% column_to_rownames(., var = "SampleName")
a_datax$bi_sex <- ifelse(a_datax$Sex == "M", 0, 1)
a_datax$bi_group <- ifelse(a_datax$Group == "Blue", 0, 1)


#checks
all(row.names(a_datax) %in% colnames(a)) #check 
a_datax <- a_datax[colnames(a),]
a <- a[,rownames(a_datax)]
all(rownames(a_datax) == colnames(a)) #check



#normalize and filter with all groups -
# I normalize with all groups since I want to compare them after
d = apply(a, 2, as.numeric)
dim(d)
# [1] 34039    36
d0= DGEList(d, group = a_datax$Sex)
dim(d0)
rownames(d0) <- rownames(a)
d0 <- calcNormFactors(d0)

cutoff <- 14
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
# [1] 10229    24

#get just females: 
females <- a_datax %>% filter(Sex == "F")

females$bi_group -> group.dl

puberty.dl <- (females$Puberty-min(females$Puberty))/(max(females$Puberty)-min(females$Puberty))

dge.f <- dge.dl[, dge.dl$samples$group %in% c("F")]
dge.f$samples$group <- droplevels(dge.f$samples$group)
dge.f$samples$group
dge.dl<- dge.f
dge.dl$samples$group

all(row.names(females) %in% colnames(dge.dl$counts)) #check 


design.dl <- model.matrix(~group.dl+puberty.dl)
colnames(design.dl) -> mycolnames

v.dl = voom(dge.dl, design.dl, plot = T)
vfit.dl = lmFit(v.dl, design.dl)
efit.dl = eBayes(vfit.dl)



p.dl.limma = efit.dl[["p.value"]]


saveRDS(v.dl, "females/clean_data/limma_vdl_females_ARC_puberty.RDS")


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

saveRDS(q.dl,("females/clean_data/limma_vdl_eFDRcutoffs_females_ARC_puberty.RDS"))


tmp1 <- contrasts.fit(efit.dl, coef = 2)#group
tmp2 <- contrasts.fit(efit.dl, coef = 3)#puberty 

limma_list <- list()
topTable(tmp1, sort.by = "P", n = Inf) %>% 
  rownames_to_column('genes') -> limma_list$group

topTable(tmp2, sort.by = "P", n = Inf) %>% 
  rownames_to_column('genes') -> limma_list$puberty


saveRDS(limma_list,"females/clean_data/limma_females_ARC_puberty.RDS")


my_logFC_threshold = 0.5 ## what should this be?

limma_list<- readRDS("females/clean_data/limma_females_ARC_puberty.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) 

#colnames(limma_list$status)
limma_list %>% 
  map(~summarise(.,Up = sum(logFC>0),
                 Down = sum(logFC<0))) %>% 
  map(~mutate(.,Total = Up + Down))




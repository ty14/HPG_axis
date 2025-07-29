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
# [1] 34039    24
d0= DGEList(d, group = a_datax$Sex)
dim(d0)
rownames(d0) <- rownames(a)
d0 <- calcNormFactors(d0)

cutoff <- 20
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
# [1] 8923   24

#get just females: 
females <- a_datax %>% filter(Sex == "F")

females$Group -> group.dl

puberty.dl <- (females$Puberty-min(females$Puberty))/(max(females$Puberty)-min(females$Puberty))

dge.f <- dge.dl[, dge.dl$samples$group %in% c("F")]
dge.f$samples$group <- droplevels(dge.f$samples$group)
dge.f$samples$group
dge.dl<- dge.f
dge.dl$samples$group

all(row.names(females) %in% colnames(dge.dl$counts)) #check 


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

saveRDS(v.dl, "females/clean_data/limma_vdl_ARC.RDS")




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

q.dl <- readRDS("females/clean_data/limma_vdl_eFDRcutoffs_ARC.RDS")
efit.dl2[["p.value"]] <- q.dl
row.names(q.dl) <- NULL
sum(duplicated(row.names(efit.dl2$coefficients)))

saveRDS(q.dl,("females/clean_data/limma_vdl_eFDRcutoffs_ARC.RDS"))

##### Analysis pulling genes out for each contrast 
tmp1 <- contrasts.fit(efit.dl2, coef = 1) #Blue vs. Orange


topTable(tmp1, sort.by = "P", n = Inf) %>% 
  rownames_to_column('genes') -> bvo


saveRDS(bvo,"females/clean_data/limma_eFDR_ARC.RDS")

fe_bvo <- readRDS("females/clean_data/limma_eFDR_ARC.RDS") 


#quick look at numbers 
fe_bvo <- readRDS("females/clean_data/limma_eFDR_ARC.RDS") %>% 
  distinct(.) %>% 
  filter(.,abs(logFC) >= 0.5) %>%
  filter(.,P.Value <0.05)

#quick look at number of genes
fe_bvo %>% filter(., P.Value<0.05) %>% 
  summarise(.,Up = sum(logFC>0.5),
            Down = sum(logFC<0.5)) %>% 
  mutate(.,Total = Up + Down)

# Up Down Total
#   Up Down Total
# 1  4    4     8

hist(fe_bvo$logFC)

##top genes 
fe_bvo_up <- fe_bvo %>% 
  filter(logFC > 0.5)

fe_bvo_up %>% arrange(-logFC) %>% head(., 20)

fe_bvo_down <- fe_bvo %>% 
  filter(logFC < 0.5)

fe_bvo_down %>% arrange(logFC) %>% head(., 20)


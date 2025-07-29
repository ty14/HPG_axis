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
head(idx)
#counts 
p <- readRDS("brain/clean_data/pit_counts.RDS") %>%
  select(-region) %>% 
  column_to_rownames(., "ensgene")


#coldata vmn
p_data <- idx
p_datax <- p_data %>% column_to_rownames(., var = "SampleName")
p_datax$bi_sex <- ifelse(p_datax$Sex == "M", 0, 1)
p_datax$bi_group <- ifelse(p_datax$Group == "Blue", 0, 1)


#checks
all(row.names(p_datax) %in% colnames(p)) #check 
p_datax <- p_datax[colnames(p),]
p <- p[,rownames(p_datax)]
all(rownames(p_datax) == colnames(p)) #check




#normalize and filter with all groups -
# I normalize with all groups since I want to compare them after
d = apply(p, 2, as.numeric)
dim(d)
# [1] 34039    47
d0= DGEList(d, group = p_datax$Sex)
dim(d0)
rownames(d0) <- rownames(p)
d0 <- calcNormFactors(d0)

cutoff <- 10
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
# [1] 10377    47

#get just males: 
males <- p_datax %>% filter(Sex == "M")

males$Group -> group.dl

puberty.dl <- (males$Puberty-min(males$Puberty))/(max(males$Puberty)-min(males$Puberty))

dge.m <- dge.dl[, dge.dl$samples$group %in% c("M")]
dge.m$samples$group <- droplevels(dge.m$samples$group)
dge.m$samples$group
dge.dl<- dge.m
dge.dl$samples$group

all(row.names(males) %in% colnames(dge.dl$counts)) #check 


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

saveRDS(v.dl, "males/clean_data/limma_vdl_PIT.RDS")




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
saveRDS(q.dl,("males/clean_data/limma_vdl_eFDRcutoffs_PIT.RDS"))

q.dl <- readRDS("males/clean_data/limma_vdl_eFDRcutoffs_PIT.RDS")
efit.dl2[["p.value"]] <- q.dl
row.names(q.dl) <- NULL
sum(duplicated(row.names(efit.dl2$coefficients)))

##### Analysis pulling genes out for each contrast 
tmp1 <- contrasts.fit(efit.dl2, coef = 1) #Blue vs. Orange


topTable(tmp1, sort.by = "P", n = Inf) %>% 
  rownames_to_column('genes') -> bvo


saveRDS(bvo,"males/clean_data/limma_eFDR_PIT.RDS")

fe_bvo <- readRDS("males/clean_data/limma_eFDR_PIT.RDS") 


#quick look at numbers 
m_bvo <- readRDS("males/clean_data/limma_eFDR_PIT.RDS") %>% 
  distinct(.) %>% 
  filter(.,abs(logFC) >= 0.5) %>%
  filter(.,P.Value <0.05)

#quick look at number of genes
m_bvo %>% filter(., P.Value<0.05) %>% 
  summarise(.,Up = sum(logFC>0.5),
            Down = sum(logFC<0.5)) %>% 
  mutate(.,Total = Up + Down)


# Up Down Total
# 7    4    11

hist(m_bvo$logFC)

##top genes 
m_bvo_up <- m_bvo %>% 
  filter(logFC > 0.5)

m_bvo_up %>% arrange(-logFC) %>% head(., 20)

m_bvo_down <- m_bvo %>% 
  filter(logFC < 0.5)

m_bvo_down %>% arrange(logFC) %>% head(., 20)




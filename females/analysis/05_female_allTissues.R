library(limma)
library(edgeR)
library(WGCNA)
library(tidyverse)

#Female WGCNA  - 


#counts 
a <- readRDS("brain/clean_data/arc_counts.RDS") %>%
  select(-region) %>% 
  column_to_rownames(., "ensgene")

colnames(a) <- paste("A", colnames(a), sep = "_")

v <- readRDS("brain/clean_data/vmn_counts.RDS") %>%
  select(-region) %>% 
  column_to_rownames(., "ensgene")

colnames(v) <- paste("V", colnames(v), sep = "_")

p <- readRDS("brain/clean_data/pit_counts.RDS") %>%
  select(-region) %>% 
  column_to_rownames(., "ensgene")

colnames(p) <- paste("P", colnames(p), sep = "_")


#gonads
g <- readRDS("raw_data/gonads_clean_counts.RDS")


# sample information 
id <- read_csv("raw_data/traitData.csv")
colnames(id)
idx <- id %>% select(1:4,15)

#coldata 
g_data <- idx
g_datax <- g_data %>% column_to_rownames(., var = "SampleName")
g_datax$bi_group <- ifelse(g_datax$Group == "Blue", 0, 1)

fg_datax <- g_datax %>% filter(Sex == "F")

#checks
all(row.names(fg_datax) %in% colnames(g)) #check 
fg <- g[,rownames(fg_datax)]
all(rownames(fg_datax) == colnames(fg)) #check

colnames(fg) <- paste("O", colnames(fg), sep = "_")

#males
mg_datax <- g_datax %>% filter(Sex == "M")

#checks
all(row.names(mg_datax) %in% colnames(g)) #check 
mg <- g[,rownames(mg_datax)]
all(rownames(mg_datax) == colnames(mg)) #check

colnames(mg) <- paste("T", colnames(mg), sep = "_")


#make coldata with id and region using colnames from each count data set 

ac <- colnames(a) %>% as.data.frame() %>% mutate(tissue = "ARC")
vc <- colnames(v) %>% as.data.frame() %>% mutate(tissue = "VMN")
pc <- colnames(p) %>% as.data.frame() %>% mutate(tissue = "PIT")
oc <- colnames(fg) %>% as.data.frame() %>% mutate(tissue = "OVT")
tc <- colnames(mg) %>% as.data.frame() %>% mutate(tissue = "TES")

col_data <- ac %>% rbind(vc,pc,oc,tc) 
colnames(col_data) <- c("id", "tissue")
col_data <- col_data %>% column_to_rownames(., var = "id")

c <- a %>% cbind(v,p, fg, mg)
all(row.names(col_data) %in% colnames(c)) #check 


# now do each tissue by group and sex 
x <- col_data %>% rownames_to_column(., var = "id") %>% 
  mutate(SampleName = substr(id, 3, 8)) %>% pivot_wider(names_from = "tissue", values_from = "id")

xx <- idx %>% select(SampleName, Sex, Group, Puberty) %>% full_join(x) %>% 
  pivot_longer(cols = 5:9, names_to = "tissue", values_to = "id") %>% na.omit %>% 
  column_to_rownames(., var = "id")

f_data <- xx %>% filter(Sex == "F")


all(row.names(f_data) %in% colnames(c)) #check 
fc <- c[,rownames(f_data)]
all(rownames(f_data) == colnames(fc)) #check


#normalize and filter with all groups -
# I normalize with all groups since I want to compare them after
d = apply(fc, 2, as.numeric)
dim(d)
# [1] 34039    74
d0= DGEList(d, group = f_data$Group)
dim(d0)
rownames(d0) <- rownames(fc)
d0 <- calcNormFactors(d0)

cutoff <- 25
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
# [1] 12623    74 #cutoff 10 
# [1] 9332   74 # cutoff 25

f_data$Group %>%
  factor(.,levels = c("Blue", "Orange")) -> group.dl


design.dl <- model.matrix(~ 0 + group.dl)
colnames(design.dl) -> mycolnames

v.dl = voom(dge.dl, design.dl, plot = T)


# Many functions expect the matrix to be transposed
datExpr <- t(v.dl$E)
## check rows/cols
nrow(datExpr)
ncol(datExpr)
rownames(datExpr)


# getting trait data
head(f_data)

trait <- f_data %>%
  mutate(Group = ifelse(.$Group == "Orange", 0, 1)) %>% 
  mutate(tissue = ifelse(.$tissue == "ARC", 0, .$tissue)) %>% 
  mutate(tissue = ifelse(.$tissue == "VMN", 1, .$tissue)) %>% 
  mutate(tissue = ifelse(.$tissue == "PIT", 2, .$tissue)) %>% 
  mutate(tissue = ifelse(.$tissue == "OVT", 3, .$tissue)) %>% 
  select(Group, Puberty, tissue)

#make everything numeric 
trait[1:3] <- lapply(trait[1:3], as.numeric)


sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)

par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 250, col = "red");
# Determine250uster under the line
clust = cutreeStatic(sampleTree, cutHeight = 60, minSize = 10)
table(clust)

sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(trait, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(trait),
                    main = "Female dendrogram and trait heatmap")

collectGarbage()
#not going to filter because ascenders rreally are not grouping anyways?
saveRDS(datExpr,"clean_data/WGCNA/WGCNA_datExpr_females_allTissues.RDS")



# Run WGCNA for each dom, descender, controls 
# Choose a set of soft-thresholding powers
powers = c(c(1:20), seq(from = 12, to=40, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, )
# Plot the results:
print(sft)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;


# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main ="Female: Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");


# this line corresponds to using an R^2 cut-off of h

abline(h=0.85,col="red")  #changed to 0.8 (originally 0.9 in the tutorial)
abline(h=0.90,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = "Female: Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
print(sft)

dev.off()

cor <- WGCNA::cor
WGCNA_get_net <- function(my_tissue = "Ovary",
                          my_power =  5,
                          my_TOMType ="signed", 
                          my_networkType = "signed hybrid"){
  
  x <- readRDS("clean_data/WGCNA/WGCNA_datExpr_females_allTissues.RDS") 
  set.seed(312)
  net = blockwiseModules(x, 
                         power = my_power,
                         TOMType = my_TOMType, 
                         networkType = my_networkType,
                         minModuleSize = 50,
                         reassignThreshold = 0, 
                         mergeCutHeight = 0.25,
                         numericLabels = TRUE, 
                         pamRespectsDendro = FALSE,
                         saveTOMs = FALSE,
                         verbose = 3)
  
  
  saveRDS(net, "clean_data/WGCNA/WGCNA_net_females_all_Power1_cutoff50.RDS")
  
}

WGCNA_get_net("Ovary", 1, "signed", "signed hybrid")
# ========================================================================================

net <- readRDS(glue::glue("clean_data/WGCNA/WGCNA_net_females_all_Power1_cutoff50.RDS"))

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath

dev.off() # make sure you do this before AND after 
png(file = "imgs/cluster_dendo_females_all_Power1_minmod50.png",
    width=600, height=350)
plotDendroAndColors(net$dendrograms[[2]], mergedColors[net$blockGenes[[2]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "Females: Power1")

dev.off()

MEs = net$MEs

##################################
set.seed(312)
datExpr <- readRDS("clean_data/WGCNA/WGCNA_datExpr_females_allTissues.RDS") 
net <- readRDS("clean_data/WGCNA/WGCNA_net_females_all_Power1.RDS")
dim(datExpr)


#PART 1
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
moduleNumber = length(unique(moduleColors))
rownames(datExpr)


# getting trait data
traitx <- trait 
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(traitx, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(traitx),
                    main = "")

datTraits <-traitx

#PART 2
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

saveRDS(MEs,"clean_data/WGCNA/WGCNA_MEs_females_all_Power1.RDS")

#PART3
# sizeGrWindow(10,6)

dev.off() # make sure you do this before AND after

png(file = "imgs/module_trait_females_all_Power1.png",
    width=2000, height=2600, res = 300)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = "Females: Module-trait relationships")
#
## remember - they are correlation, not regression

dev.off() 

#PART 4 =====================================================================================
# Define variable David's score containing the David's score column of datTrait
status = as.data.frame(datTraits$Group);
names(status) = "trait"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, status, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(status), sep="");
names(GSPvalue) = paste("p.GS.", names(status), sep="");


geneModuleMembership %>% 
  rownames_to_column("ensgene") %>% 
  left_join(geneTraitSignificance %>% rownames_to_column("ensgene")) -> gene_MM_TS_all



# to make holistic freaking dataframe ==================================================
gene_MM_TS_all$module <- moduleColors


get_mm <- function(x){
  x$moduleMembership <- x[colnames(x) == paste("MM",x$module,sep = "")] %>%
    unlist %>%
    as.numeric
  xx <- x %>%
    dplyr::select(ensgene,module,moduleMembership,GS.trait)
  return(x)
}


wgcna_whole <- get_mm(gene_MM_TS_all[1,])

for(i in 2:nrow(gene_MM_TS_all)){
  wgcna_whole <- rbind(wgcna_whole,get_mm(gene_MM_TS_all[i,]))
}

wgcna_whole %>%
  rename(GS.status = GS.trait) -> wgcna_whole


wgcna_whole %>%
  left_join(geneTraitSignificance %>% rownames_to_column("ensgene")) -> wgcna_all

saveRDS(wgcna_all, "clean_data/WGCNA/WGCNA_WGCNA_MM_GS_females_all_Power1.RDS")



##boxplots

traitxx <- xx %>% filter(Sex == "F") %>% rownames_to_column(., var = "SampleID")

ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = 2:8, names_to = "Module") %>% 
  full_join(traitxx) %>% filter(Module != "MEgrey") 

head(ME_df)

ME_df$Module <- gsub("ME", "", ME_df$Module)
source("functions/geom_boxjitter.R")


p1 <- ME_df %>% 
  ggplot(aes(tissue, value, fill = factor(interaction(tissue,Group)),color = factor(interaction(tissue,Group))))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(8))+
  scale_fill_manual(values = viridis::viridis(8))+        
  labs(x = "",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol = 6) + theme_classic()+
  theme(text =element_text(size =10))

p1

ggsave("imgs/boxplots_WGCNA_females_alltissues.png", p1, width= 15, height = 5, dpi = 150)

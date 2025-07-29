library(WGCNA)
library(limma)
library(edgeR)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)


# Expression values

df <- read.table("raw_data/all_gonad_counts_Star.gff", col.names = c("symbol", "count","id"))
head(df)
dfx <- df %>% pivot_wider(names_from = id, values_from = count) %>% column_to_rownames(., var = "symbol")
colnames(dfx) <- substr(colnames(dfx), 7, 12)
# dfxx <- dfx %>% rownames_to_column(., var = "symbol")
# saveRDS(dfxx, "raw_data/gonads_clean_counts.RDS")

#counts
dlNorm <- dfx

# trait data 
data <- read.csv("females/clean_data/female_trait_data.csv")
data <- na.omit(data) %>% column_to_rownames(., var = "SampleName") %>% 
  dplyr::select(-Sex, -Group)

data$treatment <- ifelse(data$treatment == "Control", 0, 1)

all(row.names(data) %in% colnames(dlNorm)) #check 
dlNorm <- dlNorm[,rownames(data)]
all(rownames(data) == colnames(dlNorm)) #



# remove zeros
dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]

d = apply(dlNorm, 2, as.numeric)
dim(d)
# 24287    23
d0= DGEList(d, group = data$treatment)
dim(d0)
rownames(d0) <- rownames(dlNorm)
d0 <- calcNormFactors(d0)

cutoff <- 10
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
# [1] 10270    23


data$treatment %>% as.factor(.) -> group.dl


design.dl <- model.matrix(~ 0 + group.dl)
colnames(design.dl) -> mycolnames

v.dl = voom(dge.dl, design.dl, plot = T)


# Many functions expect the matrix to be transposed
datExpr <- t(v.dl$E)
## check rows/cols
nrow(datExpr)
ncol(datExpr)
rownames(datExpr)

#make everything numeric in trait data 
data[1:4] <- lapply(data[1:4], as.numeric)
data


sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)

par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 60, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 100, minSize = 10)
table(clust)
collectGarbage()
saveRDS(datExpr,"females/clean_data/WGCNA/WGCNA_datExpr_OVT.RDS")



# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
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
     main ="Ovaries: Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");


# this line corresponds to using an R^2 cut-off of h

abline(h=0.85,col="red")  #changed to 0.8 (originally 0.9 in the tutorial)
abline(h=0.90,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = "Ovaries: Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
print(sft)
# $powerEstimate
# [1] 4


# TOMType = "signed", networkType = "signed",
# TOMType = "unsigned", networkType = "unsigned",
# TOMType = c("signed", networkType = "signed hybrid")
cor <- WGCNA::cor
WGCNA_get_net <- function(my_tissue = "OVT",
                          my_power =  4,
                          my_TOMType ="signed", 
                          my_networkType = "signed hybrid"){
  
  x <- readRDS("females/clean_data/WGCNA/WGCNA_datExpr_OVT.RDS") 
  set.seed(312)
  net = blockwiseModules(x, 
                         power = my_power,
                         TOMType = my_TOMType, 
                         networkType = my_networkType,
                         minModuleSize = 100,
                         reassignThreshold = 0, 
                         mergeCutHeight = 0.25,
                         numericLabels = TRUE, 
                         pamRespectsDendro = FALSE,
                         saveTOMs = FALSE,
                         verbose = 3)
  
  
  saveRDS(net, "females/clean_data/WGCNA/WGCNA_net_OVT_Power4_cutoff100.RDS")
  
}

WGCNA_get_net("OVT", 4, "signed", "signed hybrid")
# ========================================================================================

net <- readRDS("females/clean_data/WGCNA/WGCNA_net_OVT_Power4_cutoff100.RDS")

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath

dev.off() # make sure you do this before AND after 
png(file = "imgs/cluster_dendo_OVT_Power4_minmod100.png",
    width=600, height=350)
plotDendroAndColors(net$dendrograms[[2]], mergedColors[net$blockGenes[[2]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "Ovaries: Power4")

dev.off()

MEs = net$MEs

##################################
set.seed(312)
datExpr <- readRDS("females/clean_data/WGCNA/WGCNA_datExpr_OVT.RDS") 
net <- readRDS("females/clean_data/WGCNA/WGCNA_net_OVT_Power4_cutoff100.RDS")
dim(datExpr)
# [1]    23 10270

#PART 1
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
moduleNumber = length(unique(moduleColors))
rownames(datExpr)


# getting trait data
traitxx <- data
traitxx$`E2 pg/ml` <- traitxx$con_pg_ml
traitxx <- traitxx %>% select(-con_pg_ml, -Puberty, -Cohort)
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(traitxx, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(traitxx),
                    main = "")

datTraits <-traitxx

#PART 2
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

saveRDS(MEs,"females/clean_data/WGCNA/WGCNA_MEs_OVT_Power4.RDS")


#PART3
# sizeGrWindow(10,6)

dev.off() # make sure you do this before AND after

png(file = "imgs/module_trait_OVT_females_Power4_use.png",
    width=2000, height=2600, res = 300)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(7, 9.5, 4, 4));
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
               main = "Ovaries: Module-trait relationships")
#
## remember - they are correlation, not regression

dev.off() 
#PART 4 =====================================================================================
# Define variable David's score containing the David's score column of datTrait
status = as.data.frame(datTraits$treatment);
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
  dplyr::rename(GS.status = GS.trait) -> wgcna_whole


wgcna_whole %>%
  left_join(geneTraitSignificance %>% rownames_to_column("ensgene")) -> wgcna_all

saveRDS(wgcna_all, "females/clean_data/WGCNA/WGCNA_WGCNA_MM_GS_OVT_Power4.RDS")


#plots 
df <- read.csv("females/clean_data/female_trait_data.csv")
head(df)

MEs <- readRDS("females/clean_data/WGCNA/WGCNA_MEs_OVT_Power4.RDS")

ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleName") %>%
  pivot_longer(cols = 2:26, names_to = "Module") %>% 
  full_join(df) %>% filter(Module != "MEgrey") 

head(ME_df)

ME_df$Module <- gsub("ME", "", ME_df$Module)
source("functions/geom_boxjitter.R")


p1 <- ME_df %>% 
  ggplot(aes(treatment, value, fill = Group, color = Group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("blue", "darkorange"))+
  scale_fill_manual(values = c("blue", "darkorange"))+         
  labs(x = "",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol = 5) + theme_classic()+
  theme(legend.position = "none", text =element_text(size =10))

p1



bl <- ME_df %>% filter(Module == "black") %>% 
  ggplot(aes(treatment, value, fill = Group, color = Group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("blue", "darkorange"))+
  scale_fill_manual(values = c("blue", "darkorange"))+         
  labs(x = "",
       y = "Module eigengene")+  theme_classic()+ ggtitle("black: 350 genes")+
  theme(legend.position = "none", text =element_text(size =10))

ggsave("imgs/blackOVTmodule.png", bl, width = 2.5, height = 2.5)



blue <- ME_df %>% filter(Module == "blue") %>% 
  ggplot(aes(treatment, value, fill = Group, color = Group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("blue", "darkorange"))+
  scale_fill_manual(values = c("blue", "darkorange"))+         
  labs(x = "",
       y = "Module eigengene")+  theme_classic()+ ggtitle("blue: 952 genes")+
  theme(legend.position = "none", text =element_text(size =10))

ggsave("imgs/blueOVTmodule.png", blue, width = 2.5, height = 2.5)

yellow <- ME_df %>% filter(Module == "yellow") %>% 
  ggplot(aes(treatment, value, fill = Group, color = Group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("blue", "darkorange"))+
  scale_fill_manual(values = c("blue", "darkorange"))+         
  labs(x = "",
       y = "Module eigengene")+  theme_classic()+ ggtitle("yellow: 804 genes")+
  theme(legend.position = "none", text =element_text(size =10))

ggsave("imgs/yellowOVTmodule.png", yellow, width = 2.5, height = 2.5)

yellow


sal <- ME_df %>% filter(Module == "salmon") %>% 
  ggplot(aes(treatment, value, fill = Group, color = Group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("blue", "darkorange"))+
  scale_fill_manual(values = c("blue", "darkorange"))+         
  labs(x = "",
       y = "Module eigengene")+  theme_classic()+ ggtitle("salmon: 237 genes")+
  theme(legend.position = "none", text =element_text(size =10))

ggsave("imgs/salmonOVTmodule.png", sal, width = 2.5, height = 2.5)

sal


t <- ME_df %>% filter(Module == "brown") %>% 
  ggplot(aes(treatment, value, fill = Group, color = Group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("blue", "darkorange"))+
  scale_fill_manual(values = c("blue", "darkorange"))+         
  labs(x = "",
       y = "Module eigengene")+  theme_classic()+ ggtitle("brown: 806 genes")+
  theme(legend.position = "none", text =element_text(size =10))

ggsave("imgs/brown2OVTmodule.png", t, width = 2.5, height = 2.5)


ME_dfx <- ME_df %>% filter(Module %in% c("blue", "yellow", "midnightblue", "salmon", "brown", "turquoise", "black", "lightgreen"))

ME_dfx$Module <- factor(ME_dfx$Module, levels = c("blue", "brown", "midnightblue", "salmon", "yellow", "black", "turquoise", "lightgreen"))

p2 <- ME_df %>% filter(Module %in% c("blue", "yellow", "midnightblue", "salmon", "brown", "turquoise", "darkgreen", "black", "lightgreen")) %>% 
  ggplot(aes(treatment, value, fill = Group, color = Group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("blue", "darkorange"))+
  scale_fill_manual(values = c("blue", "darkorange"))+         
  labs(x = "",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol = 3) + theme_classic()+
  theme(legend.position = "none", text =element_text(size =10))

p2

p22 <- ME_dfx %>%  
  ggplot(aes(treatment, value, fill = Group, color = Group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("blue", "darkorange"))+
  scale_fill_manual(values = c("blue", "darkorange"))+         
  labs(x = "",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol = 5) + theme_classic()+
  theme(legend.position = "none", text =element_text(size =10))

p22

ME_dfxx <- ME_df %>% filter(Module %in% c("brown", "lightgreen", "black", "darkgreen", "cyan"))
ME_dfx$Module <- factor(ME_dfx$Module, levels = c("brown","cyan","darkgreen","black","lightgreen" )

p3 <- ME_dfxx %>% 
  ggplot(aes(con_pg_ml, value))+
  geom_point()+
  geom_smooth(method = "lm", se =F)+
  labs(x = "Estradiol pg/ml",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol = 5) + theme_classic()+
  theme(text =element_text(size =15))

p3

p3 <- ME_df %>% filter(Module %in% c("blue", "cyan", "brown", "turquoise", "black", "lightgreen")) %>% 
  ggplot(aes(con_pg_ml, value))+
  geom_point()+
  geom_smooth(method = "lm", se =F)+
  labs(x = "Estradiol pg/ml",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol = 3) + theme_classic()+
  theme(text =element_text(size =15))

p3

p4 <- ME_df %>% filter(Module %in% c("pink")) %>% 
  ggplot(aes(Puberty, value))+
  geom_point()+
  geom_smooth(method = "lm", se =F)+
  labs(x = "Puberty Onset",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol = 3) + theme_classic()+
  theme(text =element_text(size =15))

p4


b2 <- ME_df %>% filter(Module == "black") %>% 
  ggplot(aes(con_pg_ml, value))+
  geom_point()+
  geom_smooth(method = "lm", se =F)+
  labs(x = "Estradiol pg/ml",
       y = "Module eigengene")+
  # facet_wrap(~Module, ncol = 5) +
  theme_classic()+
  ggtitle("black: 350 genes")+
  theme(text =element_text(size =10))
b2
ggsave("imgs/blackE2OVTmodule.png", b2, width = 2.5, height = 2.5)


dg <- ME_df %>% filter(Module == "darkgreen") %>% 
  ggplot(aes(con_pg_ml, value))+
  geom_point()+
  geom_smooth(method = "lm", se =F)+
  labs(x = "Estradiol pg/ml",
       y = "Module eigengene")+
  # facet_wrap(~Module, ncol = 5) +
  theme_classic()+
  ggtitle("darkgreen: 160 genes")+
  theme(text =element_text(size =10))

ggsave("imgs/dgE2OVTmodule.png", dg, width = 2.5, height = 2.5)


c <- ME_df %>% filter(Module == "cyan") %>% 
  ggplot(aes(con_pg_ml, value))+
  geom_point()+
  geom_smooth(method = "lm", se =F)+
  labs(x = "Estradiol pg/ml",
       y = "Module eigengene")+
  # facet_wrap(~Module, ncol = 5) +
  theme_classic()+
  ggtitle("cyan: 233 genes")+
  theme(text =element_text(size =10))

ggsave("imgs/cyanE2OVTmodule.png", c, width = 2.5, height = 2.5)

bro <- ME_df %>% filter(Module == "brown") %>% 
  ggplot(aes(con_pg_ml, value))+
  geom_point()+
  geom_smooth(method = "lm", se =F)+
  labs(x = "Estradiol pg/ml",
       y = "Module eigengene")+
  # facet_wrap(~Module, ncol = 5) +
  theme_classic()+
  ggtitle("brown: 806 genes")+
  theme(text =element_text(size =10))

ggsave("imgs/brownE2OVTmodule.png",bro, width = 2.5, height = 2.5)


#PART 5 =====================================================================================

# Module of interest via significance 
my_trait = "Status"
module = "lightgreen"
module = "blue"
module = 'turquoise'
module = 'salmon'
module = 'yellow'
module = 'brown'
module = 'midnightblue'
module = 'darkgreen'
module = 'black'
module_list = c("green", "brown",  "blue", "turquoise","salmon", "yellow", "midnightblue", "darkgreen", "black")
hub_gene_list = vector('list', length = length(module_list))
names(hub_gene_list) <- module_list


column = match(module, modNames);
moduleGenes = moduleColors==module;
colnames(datExpr)[moduleColors==module] -> module_gene
# 
# 
# sizeGrWindow(7, 7);
# par(mfrow = c(1,1));
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                    abs(geneTraitSignificance[moduleGenes, 1]),
#                    xlab = paste("Module Membership in", module, "module"),
#                    ylab = glue::glue("Gene significance for {my_trait}"),
#                    main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
# datExpr[moduleGenes,]

my_MM_threshold = 0.8
my_GS_threshold = 0.2


for (module in module_list){
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  
  colnames(datExpr)[moduleColors==module] -> module_gene
  
  
  # grcm38 %>% 
  #   dplyr::select(ensgene, symbol, chr, description) %>% 
  #   filter(ensgene %in% module_gene) -> module_gene_info
  
  
  
  geneModuleMembership %>% 
    rownames_to_column("gene") %>% 
    left_join(geneTraitSignificance %>% rownames_to_column("gene")) -> gene_MM_TS
  
  
  
  gene_MM_TS %>% 
    filter(gene %in% module_gene) %>% 
    dplyr::select(gene, glue::glue("MM{module}"), GS.trait) %>% 
    filter(abs(GS.trait) >= my_GS_threshold) -> x
  
  x[x[,glue::glue("MM{module}")]>my_MM_threshold,] -> hub_genes
  
  hub_genes %>% 
    # left_join(module_gene_info) %>% 
    mutate(moduleName = glue::glue("{module}")) %>% 
    rename(moduleMembership = glue::glue("MM{module}")) -> hub_genes
  hub_gene_list[[module]] <- hub_genes
  
}



hub_gene_list %>% 
  do.call(rbind,.)%>% 
  unique() -> hubgenes_df


write.csv(hubgenes_df, "females/clean_data/WGCNA/WCGNA_hubgene_list_OVT_Power4.csv", row.names = F)


ADJ=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ, moduleColors)

Alldegrees1 %>%
  rownames_to_column('gene')-> kIN_df

xx <- cbind(gene = colnames(datExpr), module = moduleColors) %>%
  as.tibble()


kIN_df %>%
  left_join(xx) -> kIN_df

saveRDS(kIN_df,"clean_data/kIN_dataframe_OVT_power4.RDS")


set.seed(312)
library(clusterProfiler)
library(enrichplot)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

moduleNumber = length(unique(moduleColors))

modNames = substring(names(MEs), 3)

library(AnnotationDbi)
library(org.Rn.eg.db)


g2 <-dlNorm %>%   mutate(entrez = mapIds(org.Rn.eg.db, keys = rownames(dlNorm),
                                         column = "ENTREZID", keytype = "SYMBOL")) %>% 
  filter(!is.na(entrez)) %>% dplyr::select(entrez) %>% rownames_to_column(., var = "gene") 

gettop10GO_WGCNA <- function(module,my_showCategory = 10){
  
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  colnames(datExpr)[moduleColors==module] -> module_gene
  
  
  g2%>% 
    filter(gene %in% module_gene) %>% 
    filter(!is.na(entrez)) %>% 
    dplyr::select(entrez) -> go_df_wgcna
  
  
  ggo <- enrichGO(gene = go_df_wgcna$entrez %>% unique(),
                  OrgDb = org.Rn.eg.db::org.Rn.eg.db,
                  keyType = "ENTREZID",
                  ont = 'BP',
                  readable = T,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.50)
  
  
  fortify(
    ggo,
    showCategory = my_showCategory,
    by = "Count",
    split = NULL,
    includeAll = TRUE
  ) %>% 
    dplyr::arrange(desc(GeneRatio)) %>% 
    mutate(module = module) -> temp1
  
  return(rbind(temp1))
  
}

my_ont = "BP"
my_showCategory = 100


module_list %>% unique() -> allcolors

WGCNA_GOs <- vector('list', length(allcolors))

for(i in 1:length(allcolors)){
  gettop10GO_WGCNA(allcolors[i],my_showCategory) -> WGCNA_GOs[[i]]
}

WGCNA_GOs %>% 
  do.call(rbind,.) -> wgcna_all_gos



write.csv(wgcna_all_gos, 
          "females/clean_data/WGCNA/WGCNA_all_gos_catogeryBP_OVT_Power4.csv",
          row.names = F)

##### linear models 

ME_df

blue <- ME_df %>% filter(Module == "blue")
bl <- lmerTest::lmer(value ~treatment+con_pg_ml +(1|Cohort), data = blue)
summary(bl)


yellow <- ME_df %>% filter(Module == "yellow")
yl <- lmerTest::lmer(value ~treatment+con_pg_ml +(1|Cohort), data = yellow)
summary(yl)

# Brown is E2 and group
brown <- ME_df %>% filter(Module == "brown")
brl <- lmerTest::lmer(value ~treatment+con_pg_ml +(1|Cohort), data = brown)
summary(brl)

tu <- ME_df %>% filter(Module == "turquoise")
tul <- lmerTest::lmer(value ~treatment+con_pg_ml +(1|Cohort), data = tu)
summary(tul)


mid <- ME_df %>% filter(Module == "midnightblue")
mbl <- lmerTest::lmer(value ~treatment+con_pg_ml+(1|Cohort), data = mid)
summary(mbl)

s <- ME_df %>% filter(Module == "salmon")
sl <- lmerTest::lmer(value ~treatment+con_pg_ml+(1|Cohort), data = s)
summary(sl)

# Black is E2 and group
b <- ME_df %>% filter(Module == "black")
bl <- lmerTest::lmer(value ~treatment+con_pg_ml +(1|Cohort), data = b)
summary(bl)
cor.test(c$value, c$con_pg_ml)
anova(bl)
library(MuMIn)
r.squaredGLMM(bl) 

#only E2
dg <- ME_df %>% filter(Module == "darkgreen")
dgl <- lmerTest::lmer(value ~treatment+con_pg_ml+(1|Cohort), data = dg)
summary(dgl)

# Black is E2 and group
lg <- ME_df %>% filter(Module == "lightgreen")
lgl <- lmerTest::lmer(value ~treatment+con_pg_ml+(1|Cohort), data = lg)
summary(lgl)

# only E2
c <- ME_df %>% filter(Module == "cyan")
cl <- lmerTest::lmer(value ~treatment+con_pg_ml+(1|Cohort), data = c)
summary(cl)


### hubgenes
library(tidyverse)
my_logFC_threshold = 0.2

# females ovt
fe_ovt <- readRDS("females/clean_data/limma_eFDR_OVT.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val) %>% distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold) %>%
  filter(.,P.Value <0.05)

o_up <- fe_ovt %>% filter(logFC >= 0.2) %>% arrange(-logFC) #762
o_down <- fe_ovt %>% filter(logFC <= -0.2)%>% arrange(logFC) #533


hub <- read_csv("females/clean_data/WGCNA/WCGNA_hubgene_list_OVT_Power4.csv")
head(hub)


hb <- hub %>% filter(moduleName == "blue") %>% arrange(-moduleMembership) #157 out of 952
head(hb, 10)
hb$gene[hb$gene %in% o_down$genes] #40
hb$gene[hb$gene %in% o_up$genes] #0
# goterm
# ribosome biogenesis
# ncRNA processing
# protein-RNA complex assembly
# protein-RNA complex organization
# non-membrane-bounded organelle assembly


hbr <- hub %>% filter(moduleName == "brown") %>% arrange(-moduleMembership) #188 out of 806
head(hbr, 10)
hbr$gene[hbr$gene %in% o_down$genes] #14
hbr$gene[hbr$gene %in% o_up$genes]
#goterms
# generation of precursor metabolites and energy
# energy derivation by oxidation of organic compounds
# cellular respiration
# aerobic respiration


hy <- hub %>% filter(moduleName == "yellow") %>% arrange(-moduleMembership) #142 out of 804
head(hy, 10)

hy$gene[hy$gene %in% o_down$genes] #108
hy$gene[hy$gene %in% o_up$genes]
# goterms
# muscle cell differentiation
# muscle system process
# actin filament organization
# positive regulation of cell adhesion
# striated muscle cell differentiation

hmb <- hub %>% filter(moduleName == "midnightblue") %>% arrange(-moduleMembership) #21 out of 228
head(hmb, 10)

hmb$gene[hmb$gene %in% o_down$genes] #12
hmb$gene[hmb$gene %in% o_up$genes]
# Goterms
# leukocyte mediated immunity
# lymphocyte mediated immunity
# adaptive immune response
# cytokine-mediated signaling pathway


hs <- hub %>% filter(moduleName == "salmon") %>% arrange(-moduleMembership) #35 out of 237
head(hs, 10)
hs$gene[hs$gene %in% o_down$genes] #11
hs$gene[hs$gene %in% o_up$genes]
# goterms
# ncRNA processing
# mRNA processing
# rRNA metabolic process
# rRNA processing


ht <- hub %>% filter(moduleName == "turquoise") %>% arrange(-moduleMembership) #292 out of 1601
top <- head(ht, 10)

ht$gene[ht$gene %in% o_down$genes] #0
ht$gene[ht$gene %in% o_up$genes] #221
o_up[o_up$genes %in% top$gene,]

# goterms
# cilium assembly
# cilium organization


hbl <- hub %>% filter(moduleName == "black") %>% arrange(-moduleMembership) #42 out of 350
head(hbl, 10)
hbl$gene[hbl$gene %in% o_down$genes] #0
hbl$gene[hbl$gene %in% o_up$genes] #34
# goterms
# regulation of leukocyte apoptotic process
# leukocyte apoptotic process


hlg <- hub %>% filter(moduleName == "lightgreen") %>% arrange(-moduleMembership) #0 total greens 188
head(hlg, 10)


hdg <- hub %>% filter(moduleName == "darkgreen") %>% arrange(-moduleMembership) #5 out of 160
head(hdg, 10)

hdg$gene[hdg$gene %in% o_down$genes] #0
hdg$gene[hdg$gene %in% o_up$genes] #2
# goterms
# cellular response to peptide hormone stimulus
# response to insulin
# response to estradiol
# placenta development



hc <- hub %>% filter(moduleName == "cyan") %>% arrange(-moduleMembership) #0 total greens 233
head(hc, 10)



###module number
datExpr <- readRDS("females/clean_data/WGCNA/WGCNA_datExpr_OVT.RDS") 
net <- readRDS("females/clean_data/WGCNA/WGCNA_net_OVT_Power4_cutoff100.RDS")

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

moduleNumber = length(unique(moduleColors))

modNames = substring(names(MEs), 3)


table(module)

moduleColors %>% 
  table() %>% 
  as.data.frame() %>% arrange(Freq)   -> modnum 
colnames(modnum) <- c("module","count")
modnum <- modnum %>% filter(module !="grey")

modnum$module <- factor(modnum$module,levels=unique(modnum$module[order(modnum$count)]))
mycol <- modnum %>% 
  .$module %>% as.character() 
# for some reason this has to be character again to be ordered properly in the figure...!! 


modnum %>%
  filter(module != "grey") %>% 
  ggplot(aes(y = module, x = count, fill = module))+
  geom_bar(stat = 'identity', color = NA)+
  geom_text(aes(label = count),hjust = -0.2, size = 4)+
  xlim(0, max(modnum$count)+75)+
  scale_fill_manual(values = mycol)+
  theme_minimal(base_size = 7)+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(margin = margin(r=-6)),
        panel.grid = element_blank())+
  labs(x = "",
       y = "",title = "Module Size")+
  theme_minimal(base_size = 15)+
  theme(legend.position = "none", 
        legend.title = element_text(size =rel(1.2)),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size =rel(1.5)),
        legend.key.size = unit(.9, 'cm'),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_blank(),
        # panel.background = element_blank(),
        # panel.grid.minor = element_line(color = 'grey50'),
        panel.grid = element_blank(),
        axis.text.y = element_text(hjust=1,vjust=0.25,size = rel(1.25))) -> temp_p
temp_p


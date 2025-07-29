library(limma)
library(edgeR)
library(WGCNA)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)


# Expression values
dlNorm <-  readRDS("brain/clean_data/arc_counts.RDS") %>%
  select(-region) %>% 
  column_to_rownames(., "ensgene")

#remove zeros
dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]

#save rownames
rownames(dlNorm) -> rownames

# sample information 
id <- read_csv("raw_data/traitData.csv")
colnames(id)
idx <- id %>% select(1:4,15)
idxx <- idx %>% column_to_rownames(., var = "SampleName")

#check to see if ids match 
dlNorm <- dlNorm[,rownames(idxx)]
idxx <- idxx[colnames(dlNorm),]
dlNorm <- dlNorm[,rownames(idxx)]
all(rownames(idxx) == colnames(dlNorm)) #check
rownames(dlNorm) <- rownames

#use all data for now. 
# #getting just female data 
# females <- idxx %>% filter(Sex == "F")
# 
# #check to see if ids match 
# dlNorm_f <- dlNorm[,rownames(females)]
# all(rownames(females) == colnames(dlNorm_f)) #check

#normalize and filter with all groups -
# I normalize with all groups since I want to compare them after
rownames(dlNorm) <- rownames

d = apply(dlNorm, 2, as.numeric)
dim(d)
# [1] 23005    24
d0= DGEList(d, group = idxx$Group)
dim(d0)
rownames(d0) <- rownames(dlNorm)
d0 <- calcNormFactors(d0)

cutoff <- 20
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
# [1] 8923   24

idxx$Group %>%
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
head(idxx)

trait <- idxx %>% 
  mutate(Sex = ifelse(.$Sex == "F", 1, 0)) %>% 
  mutate(Group = ifelse(.$Group == "Orange", 0, 1)) 

#make everything numeric 
trait[1:2] <- lapply(trait[1:2], as.numeric)


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
clust = cutreeStatic(sampleTree, cutHeight = 60, minSize = 10)
table(clust)

#outliers: NMX03B male blue group (1X)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExprx = datExpr[keepSamples, ]
nGenes = ncol(datExprx)
nSamples = nrow(datExprx)

traitx <- trait[c(1:11, 13:24),]

sampleTree2 = hclust(dist(datExprx), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(traitx, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(traitx),
                    main = "VMN Males and Females dendrogram and trait heatmap")

collectGarbage()
# saveRDS(datExprx,"clean_data/WGCNA/WGCNA_datExpr_Fe_Male_VMN.RDS")



# Run WGCNA for each dom, descender, controls 
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
     main ="VMN: Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");


# this line corresponds to using an R^2 cut-off of h

abline(h=0.85,col="red")  #changed to 0.8 (originally 0.9 in the tutorial)
abline(h=0.90,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = "VMN: Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
print(sft)
# $powerEstimate
# [1] 8


dev.off()

# TOMType = "signed", networkType = "signed",
# TOMType = "unsigned", networkType = "unsigned",
# TOMType = c("signed", networkType = "signed hybrid")
cor <- WGCNA::cor
WGCNA_get_net <- function(my_tissue = "VMN",
                          my_power =  8,
                          my_TOMType ="signed", 
                          my_networkType = "signed hybrid"){
  
  x <- readRDS("clean_data/WGCNA/WGCNA_datExpr_Fe_Male_VMN.RDS") 
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
  
  
  saveRDS(net, "clean_data/WGCNA/WGCNA_net_VMN_FM_Power8.RDS")
  
}

WGCNA_get_net("VMN", 8, "signed", "signed hybrid")
# ========================================================================================

net <- readRDS(glue::glue("clean_data/WGCNA/WGCNA_net_VMN_FM_Power8.RDS"))

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath

dev.off() # make sure you do this before AND after 
png(file = "imgs/cluster_dendo_FM_VMN_Power8_minmod50.png",
    width=600, height=350)
plotDendroAndColors(net$dendrograms[[2]], mergedColors[net$blockGenes[[2]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "Ovaries: Power5")

dev.off()

MEs = net$MEs



##################################
set.seed(312)
datExpr <- readRDS("clean_data/WGCNA/WGCNA_datExpr_Fe_Male_VMN.RDS") 
net <- readRDS("clean_data/WGCNA/WGCNA_net_VMN_FM_Power8.RDS")
dim(datExpr)


#PART 1
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
moduleNumber = length(unique(moduleColors))
rownames(datExpr)


# getting trait data
traitx <- trait[c(1:11, 13:24),]
traitxx <- traitx %>% select(Group,Sex)
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(traitxx, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(traitxx),
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

saveRDS(MEs,"clean_data/WGCNA/WGCNA_MEs_FM_VMN_Power8.RDS")

#PART3
# sizeGrWindow(10,6)

dev.off() # make sure you do this before AND after

png(file = "imgs/module_trait_VMN_Power8.png",
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
               main = "Ovaries: Module-trait relationships")
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

saveRDS(wgcna_all, "clean_data/WGCNA/WGCNA_WGCNA_MM_GS_FM_VMN_Power8.RDS")



##boxplots

tdata <-idxx %>% rownames_to_column(., var = "SampleID") %>% filter(SampleID != "NMX03B")
ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = 2:14, names_to = "Module") %>% 
  full_join(tdata) %>% filter(Module != "MEgrey") 

head(ME_df)

ME_df$Module <- gsub("ME", "", ME_df$Module)
source("functions/geom_boxjitter.R")


p1 <- ME_df %>% 
  ggplot(aes(Sex, value, fill = Group, color = Group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("blue", "darkorange"))+
  scale_fill_manual(values = c("blue", "darkorange"))+         
  labs(x = "",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol = 5) + theme_classic()+
  theme(text =element_text(size =10))

p1

# ggsave("imgs/boxplots_WGCNA_females_ALL.png", p1, width= 15, height = 15, dpi = 150)


p2 <- ME_df %>% filter(Module %in% c("brown", "turquoise", "magenta", "black", "cyan", "blue", "green")) %>% 
  ggplot(aes(Group, value, fill = Group, color = Group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("blue", "darkorange"))+
  scale_fill_manual(values = c("blue", "darkorange"))+         
  labs(x = "",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol = 4) + theme_classic()+
  theme(legend.position = "none", text =element_text(size =15))

p2


ggsave("imgs/boxplots_WGCNA_females_interesting.png", p2, width= 12, height = 6, dpi = 150)


p1 <- ME_df %>% 
  ggplot(aes(Puberty, value, color = interaction(Sex, Group)))+
  geom_point() +
  geom_smooth(method = "lm", se = F)+
  scale_color_manual(values = c("blue",  "lightblue","darkorange", "pink"))+
  scale_fill_manual(values = c("blue",  "lightblue","darkorange", "pink"))+         
  labs(x = "",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol = 5) + theme_classic()+
  theme(text =element_text(size =10))

p1

#PART 5 =====================================================================================

# Module of interest via significance 
my_trait = "Status"
module = "pink"
module = "blue"
module = 'turquoise'
module = 'red'
module = 'yellow'
module = 'brown'
module = 'magenta'

module_list = c("pink", "brown",  "blue", "turquoise","red", "yellow", "magenta")
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


write.csv(hubgenes_df, "clean_data/WGCNA/WCGNA_hubgene_list_FM_VMN_Power8.csv", row.names = F)


ADJ=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ, moduleColors) 

Alldegrees1 %>% 
  rownames_to_column('gene')-> kIN_df

xx <- cbind(gene = colnames(datExpr), module = moduleColors) %>% 
  as.tibble()


kIN_df %>% 
  left_join(xx) -> kIN_df

saveRDS(kIN_df,"clean_data/kIN_dataframe_FM_VMN_power8.RDS")


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
  filter(!is.na(entrez)) %>% select(entrez) %>% rownames_to_column(., var = "gene") 

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
          "clean_data/WGCNA/WGCNA_all_gos_catogeryBP_FM_VMN_Power8.csv",
          row.names = F)




######linear model ## come back to this later. 

traitxx <- idxx %>% filter(Sex == "F") %>% rownames_to_column(., var = "SampleID")
ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = 2:21, names_to = "Module") %>% 
  full_join(traitxx) %>% filter(Module != "MEgrey") 

head(ME_df)

ME_df$Module <- gsub("ME", "", ME_df$Module)


#brown
brown <- ME_df %>% filter(Module == "brown")

mb <- lm(value~Group, data = brown)
summary(mb)
# Group Orange  0.16717    0.08347   2.003   0.0583 .

#turquoise
tur <- ME_df %>% filter(Module == "turquoise")

mt <- lm(value~Group, data = tur)
summary(mt)
# GroupOrange  0.20416    0.07945   2.570   0.0179 *


#magenta
mag <- ME_df %>% filter(Module == "magenta")

mm <- lm(value~Group, data = mag)
summary(mm)
# GroupOrange  0.16477    0.08369   1.969   0.0623 .


#magenta
black<- ME_df %>% filter(Module == "black")

mbl <- lm(value~Group, data = black)
summary(mbl)
# GroupOrange -0.22624    0.07655  -2.955  0.00755 **

#cyan
cyan<- ME_df %>% filter(Module == "cyan")

mc <- lm(value~Group, data = cyan)
summary(mc)
# GroupOrange -0.23808    0.07482  -3.182  0.00449 **


#blue
blue<- ME_df %>% filter(Module == "blue")

mblue <- lm(value~Group, data = blue)
summary(mblue)
# GroupOrange -0.25096    0.07279  -3.448  0.00241 **


#green
g<- ME_df %>% filter(Module == "green")

mg <- lm(value~Group, data = g)
summary(mg)
# GroupOrange -0.18134    0.08204  -2.210   0.0383 *


#########module number. 


datExpr <- readRDS("clean_data/WGCNA/WGCNA_datExpr_females.RDS") 
net <- readRDS("clean_data/WGCNA/WGCNA_net_females_Power5.RDS")

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


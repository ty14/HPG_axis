library(WGCNA)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)


# Expression values
#vmn counts 
v_dlNorm <-  readRDS("brain/clean_data/vmn_counts.RDS") %>%
  dplyr::select(-region) %>% 
  column_to_rownames(., "ensgene")

colnames(v_dlNorm) <- paste("V", colnames(v_dlNorm), sep = "_")

#arc counts
a_dlNorm <-  readRDS("brain/clean_data/arc_counts.RDS") %>%
  dplyr::select(-region) %>% 
  column_to_rownames(., "ensgene")

colnames(a_dlNorm) <- paste("A", colnames(a_dlNorm), sep = "_")


#pit counts
p_dlNorm <-  readRDS("brain/clean_data/pit_counts.RDS") %>%
  dplyr::select(-region) %>% 
  column_to_rownames(., "ensgene")


colnames(p_dlNorm) <- paste("P", colnames(p_dlNorm), sep = "_")

#trait data females 
#VMN
data <- read.csv("females/clean_data/female_trait_data.csv")
data <- na.omit(data)

data$SampleName <- paste("V", data$SampleName, sep = "_")
v_data <- data %>% column_to_rownames(., var = "SampleName")


all(rownames(v_data) %in% colnames(v_dlNorm)) #check 
v_data <- v_data[colnames(v_dlNorm),]
v_data <- na.omit(v_data)
v_dlNorm <- v_dlNorm[,rownames(v_data)]
all(rownames(v_data) == colnames(v_dlNorm)) #check

#ARC
data <- read.csv("females/clean_data/female_trait_data.csv")
data <- na.omit(data)

data$SampleName <- paste("A", data$SampleName, sep = "_")
a_data <- data %>% column_to_rownames(., var = "SampleName")

all(rownames(a_data) %in% colnames(a_dlNorm)) #check 
a_data <- a_data[colnames(a_dlNorm),]
a_data <- na.omit(a_data)
a_dlNorm <- a_dlNorm[,rownames(a_data)]
all(rownames(a_data) == colnames(a_dlNorm)) #check


#PIT
data <- read.csv("females/clean_data/female_trait_data.csv")
data <- na.omit(data)

data$SampleName <- paste("P", data$SampleName, sep = "_")
p_data <- data %>% column_to_rownames(., var = "SampleName")

all(row.names(p_data) %in% colnames(p_dlNorm)) #check 
# p_data <- p_data[colnames(p_dlNorm),]
p_dlNorm <- p_dlNorm[,rownames(p_data)]
all(rownames(p_data) == colnames(p_dlNorm)) #check


#bring everything together
dlNorm <- v_dlNorm %>% cbind(a_dlNorm, p_dlNorm)

data <- v_data %>% rbind(a_data, p_data)

data$treatment <- ifelse(data$treatment == "Control", 0, 1)
datax <- data %>% dplyr::select(-Sex, -Group)

dlNorm <- dlNorm[,rownames(datax)]
all(rownames(datax) == colnames(dlNorm)) #check



# remove zeros
dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]

d = apply(dlNorm, 2, as.numeric)
dim(d)
# [1]24342    51
d0= DGEList(d, group = datax$treatment)
dim(d0)
rownames(d0) <- rownames(dlNorm)
d0 <- calcNormFactors(d0)

cutoff <- 10
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
# [1] 11722    51


datax$treatment %>% as.factor(.) -> group.dl


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
datax[1:4] <- lapply(datax[1:4], as.numeric)
datax


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
saveRDS(datExpr,"females/clean_data/WGCNA/WGCNA_datExpr_brainPit.RDS")



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
WGCNA_get_net <- function(my_tissue = "VAP",
                          my_power =  8,
                          my_TOMType ="signed", 
                          my_networkType = "signed hybrid"){
  
  x <- readRDS("females/clean_data/WGCNA/WGCNA_datExpr_brainPit.RDS") 
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
  
  
  saveRDS(net, "females/clean_data/WGCNA/WGCNA_net_brainPit_Power8.RDS")
  
}

WGCNA_get_net("VAP", 8, "signed", "signed hybrid")
# ========================================================================================

net <- readRDS("females/clean_data/WGCNA/WGCNA_net_brainPit_Power8.RDS")

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath

dev.off() # make sure you do this before AND after 
png(file = "imgs/cluster_dendo_VAP_Power8_minmod50.png",
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
datExpr <- readRDS("females/clean_data/WGCNA/WGCNA_datExpr_brainPit.RDS") 
net <- readRDS("females/clean_data/WGCNA/WGCNA_net_brainPit_Power8.RDS")
dim(datExpr)
# [1]    51 11722

#PART 1
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
moduleNumber = length(unique(moduleColors))
rownames(datExpr)


# getting trait data
traitxx <- datax
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

saveRDS(MEs,"females/clean_data/WGCNA/WGCNA_MEs_VAP_Power8.RDS")


#PART3
# sizeGrWindow(10,6)

dev.off() # make sure you do this before AND after

png(file = "imgs/module_trait_VAP_females_Power8.png",
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
  rename(GS.status = GS.trait) -> wgcna_whole


wgcna_whole %>%
  left_join(geneTraitSignificance %>% rownames_to_column("ensgene")) -> wgcna_all

saveRDS(wgcna_all, "females/clean_data/WGCNA/WGCNA_WGCNA_MM_GS_VAP_Power8.RDS")


#plots 
datay <- data %>% rownames_to_column(., "SampleName") 
datay$region <- ifelse(grepl("P", datay$SampleName), "PIT", datay$SampleName)
datay$region <- ifelse(grepl("A", datay$region), "ARC", datay$region)
datay$region <- ifelse(grepl("V", datay$region), "VMN", datay$region)
datay$treatmentx <- ifelse(datay$treatment == 0, "Control", "NMX")

ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleName") %>%
  pivot_longer(cols = 2:8, names_to = "Module") %>% 
  full_join(datay) %>% filter(Module != "MEgrey") 

head(ME_df)

ME_df$Module <- gsub("ME", "", ME_df$Module)
source("functions/geom_boxjitter.R")


p1 <- ME_df %>% 
  ggplot(aes(treatmentx, value, fill = Group, color = Group))+
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

p1

# ggsave("imgs/boxplots_WGCNA_females_ALL.png", p1, width= 15, height = 15, dpi = 150)


p2 <- ME_df %>% 
  ggplot(aes(region, value, fill = region, color = region))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("green", "blue", "red"))+
  scale_fill_manual(values = c( "green", "blue", "red"))+         
  labs(x = "",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol = 3) + theme_classic()+
  theme(legend.position = "none", text =element_text(size =15))

p2


p3 <- ME_df %>% 
  ggplot(aes(con_pg_ml, value, color = region))+
  geom_point()+
  geom_smooth(method = "lm", se =F)+
  labs(x = "Estradiol pg/ml",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol = 3) + theme_classic()+
  theme(text =element_text(size =15))

p3

p4 <- ME_df %>% 
  ggplot(aes(Puberty, value, color = region))+
  geom_point()+
  geom_smooth(method = "lm", se =F)+
  labs(x = "Puberty onset",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol = 3) + theme_classic()+
  theme(text =element_text(size =15))

p4


#PART 5 =====================================================================================

# Module of interest via significance 
my_trait = "Status"
module = "green"
module = "blue"
module = 'turquoise'
module = 'red'
module = 'yellow'
module = 'brown'
module_list = c("green", "brown",  "blue", "turquoise","red", "yellow")
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


# write.csv(hubgenes_df, "females/clean_data/WGCNA/WCGNA_hubgene_list_VAP_Power8.csv", row.names = F)


# ADJ=abs(cor(datExpr,use="p"))^6
# Alldegrees1=intramodularConnectivity(ADJ, moduleColors) 
# 
# Alldegrees1 %>% 
#   rownames_to_column('gene')-> kIN_df
# 
# xx <- cbind(gene = colnames(datExpr), module = moduleColors) %>% 
#   as.tibble()
# 
# 
# kIN_df %>% 
#   left_join(xx) -> kIN_df
# 
# saveRDS(kIN_df,"clean_data/kIN_dataframe_FM_VMN_power8.RDS")


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
          "females/clean_data/WGCNA/WGCNA_all_gos_catogeryBP_VAP_Power8.csv",
          row.names = F)



#############linear models 
head(ME_df)
ME_df$regionx <- factor(ME_df$region, levels = c("PIT", "VMN", "ARC"))

b <- ME_df %>% filter(Module == "blue")
head(b)

library(lmerTest)
bl <- lmer(value ~region + (1|Cohort), data = b)
summary(bl)

bl2 <- lmer(value ~regionx + (1|Cohort), data = b)
summary(bl2)


t <- ME_df %>% filter(Module == "turquoise")
head(t)

tl <- lmer(value ~region + (1|Cohort), data = t)
summary(tl)

tl2 <- lmer(value ~regionx + (1|Cohort), data = t)
summary(tl2)


g <- ME_df %>% filter(Module == "green")
head(g)

gl <- lmer(value ~region + (1|Cohort), data = g)
summary(gl)

gl2 <- lmer(value ~regionx + (1|Cohort), data = g)
summary(gl2)


br <- ME_df %>% filter(Module == "brown")
head(br)

brl <- lmer(value ~region + (1|Cohort), data = br)
summary(brl)

brl2 <- lmer(value ~regionx + (1|Cohort), data = br)
summary(brl2)



yr <- ME_df %>% filter(Module == "yellow")
head(yr)

yrl <- lmer(value ~region + (1|Cohort), data = yr)
summary(yrl)

yr2 <- lmer(value ~regionx + (1|Cohort), data = yr)
summary(yr2)


r <- ME_df %>% filter(Module == "red")
head(r)

rl <- lmer(value ~region + (1|Cohort), data = r)
summary(rl)

r2 <- lmer(value ~regionx + (1|Cohort), data = r)
summary(r2)


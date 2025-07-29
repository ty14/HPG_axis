library(limma)
library(edgeR)
library(WGCNA)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)



df <- read.table("raw_data/all_gonad_counts_Star.gff", col.names = c("symbol", "count","id"))
head(df)
dfx <- df %>% pivot_wider(names_from = id, values_from = count) %>% column_to_rownames(., var = "symbol")
colnames(dfx) <- substr(colnames(dfx), 7, 12)
# dfxx <- dfx %>% rownames_to_column(., var = "symbol")
# saveRDS(dfxx, "raw_data/gonads_clean_counts.RDS")

#counts
g <- dfx

# sample information 
id <- read_csv("males/clean_data/male_trait_Data_wgcna.csv")
colnames(id)
idx <- id %>% select(1:3,6)

#coldata
g_data <- idx%>% filter(SampleName != "NMX03B")
g_datax <- g_data %>% column_to_rownames(., var = "SampleName")
g_datax$bi_group <- ifelse(g_datax$treatment == "Control", 0, 1)

mg_datax <- g_datax 

#checks
all(row.names(mg_datax) %in% colnames(g)) #check 
g <- g[,rownames(mg_datax)]
all(rownames(mg_datax) == colnames(g)) #check


#remove zeros
dlNorm <- g
dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]

#save rownames
rownames(dlNorm) -> rownames
#check to see if ids match 
dlNorm <- dlNorm[,rownames(mg_datax)]
all(rownames(mg_datax) == colnames(dlNorm)) #check
rownames(dlNorm) <- rownames


#normalize and filter with all groups -
# I normalize with all groups since I want to compare them after
d = apply(dlNorm, 2, as.numeric)
dim(d)
# [1] 26336    23
d0= DGEList(d, group = mg_datax$Group)
dim(d0)
rownames(d0) <- rownames(dlNorm)
d0 <- calcNormFactors(d0)

cutoff <- 10
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
# [1] 9613   23

mg_datax$treatment %>%
  factor(.,levels = c("NMX", "Control")) -> group.dl


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
head(mg_datax)

trait <- mg_datax %>% 
  select(bi_group, Cohort, con_ng_ml)

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
abline(h = 60, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 60, minSize = 10)
table(clust)

sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(trait, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(trait),
                    main = "Testes dendrogram and trait heatmap")

collectGarbage()
#not going to filter because ascenders rreally are not grouping anyways?
saveRDS(datExpr,"males/clean_data/WGCNA/WGCNA_datExpr_males.RDS")



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
     main ="Testes: Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");


# this line corresponds to using an R^2 cut-off of h

abline(h=0.85,col="red")  #changed to 0.8 (originally 0.9 in the tutorial)
abline(h=0.90,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = "Testes: Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
print(sft)
# $powerEstimate
# [1] 4


dev.off()

cor <- WGCNA::cor
WGCNA_get_net <- function(my_tissue = "Testes",
                          my_power =  4,
                          my_TOMType ="signed", 
                          my_networkType = "signed hybrid"){
  
  x <- readRDS("males/clean_data/WGCNA/WGCNA_datExpr_males.RDS") 
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
  
  
  saveRDS(net, "males/clean_data/WGCNA/WGCNA_net_males_Power4.RDS")
  
}

WGCNA_get_net("Testes", 4, "signed", "signed hybrid")


net <- readRDS(glue::glue("males/clean_data/WGCNA/WGCNA_net_males_Power4.RDS"))

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath

dev.off() # make sure you do this before AND after 
png(file = "imgs/cluster_dendo_males_Power4_minmod100.png",
    width=600, height=350)
plotDendroAndColors(net$dendrograms[[2]], mergedColors[net$blockGenes[[2]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "Testes: Power4")

dev.off()

MEs = net$MEs


##################################
set.seed(312)
datExpr <- readRDS("males/clean_data/WGCNA/WGCNA_datExpr_males.RDS") 
net <- readRDS("males/clean_data/WGCNA/WGCNA_net_males_Power4.RDS")
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

datTraits <-traitx %>% select(-Cohort, treatment=bi_group,con_ng_ml)

#PART 2
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

saveRDS(MEs,"males/clean_data/WGCNA/WGCNA_MEs_males_Power4.RDS")

#PART3
# sizeGrWindow(10,6)

dev.off() # make sure you do this before AND after

png(file = "imgs/module_trait_males_Power4.png",
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
               main = "Testes: Module-trait relationships")
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

saveRDS(wgcna_all, "males/clean_data/WGCNA/WGCNA_WGCNA_MM_GS_males_Power4.RDS")




traitx <- trait %>% rownames_to_column(., var ="SampleName")
ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleName") %>%
  pivot_longer(cols = 2:11, names_to = "Module") %>% 
  full_join(traitx) %>% filter(Module != "MEgrey") 

head(ME_df)

ME_df$Module <- gsub("ME", "", ME_df$Module)
source("functions/geom_boxjitter.R")

ME_df$treatment <- ifelse(ME_df$bi_group == 0, "Control", "NMX")

p1 <- bi_groupp1 <- ME_df %>% 
  ggplot(aes(treatment, value, fill = treatment, color = treatment))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c( "darkorange", "blue"))+
  scale_fill_manual(values = c("darkorange", "blue"))+         
  labs(x = "",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol = 5) + theme_classic()+
  theme(legend.position = "none", text =element_text(size =10))

p1

ggsave("imgs/boxplots_WGCNA_males_ALL.png", p1, width= 15, height = 15, dpi = 150)


p1 <- ME_df %>% 
  ggplot(aes(con_ng_ml, value, fill = treatment, color = treatment))+
  geom_point()+
  geom_smooth(method = "lm", se = F)+
  scale_color_manual(values = c( "darkorange", "blue"))+
  scale_fill_manual(values = c("darkorange", "blue"))+          
  labs(x = "",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol = 5) + theme_classic()+
  theme(legend.position = "none", text =element_text(size =10))

p1



t <- ME_df %>% filter(Module == "turquoise") %>% 
  ggplot(aes(treatment, value, fill = treatment, color = treatment))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c( "darkorange", "blue"))+
  scale_fill_manual(values = c("darkorange", "blue"))+          
  labs(x = "",
       y = "Module eigengene")+  theme_classic()+ ggtitle("turquoise: 1301 genes")+
  theme(legend.position = "none", text =element_text(size =10))

ggsave("imgs/tur_male_module.png", t, width = 2.5, height = 2.5)


t2 <- ME_df %>% filter(Module == "brown") %>% 
  ggplot(aes(treatment, value, fill = treatment, color = treatment))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c( "darkorange", "blue"))+
  scale_fill_manual(values = c("darkorange", "blue"))+            
  labs(x = "",
       y = "Module eigengene")+  theme_classic()+ ggtitle("brown: 840 genes")+
  theme(legend.position = "none", text =element_text(size =10))

ggsave("imgs/brown_male_module.png", t2, width = 2.5, height = 2.5)

t2



t3 <- ME_df %>% filter(Module == "blue") %>% 
  ggplot(aes(treatment, value, fill = treatment, color = treatment))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c( "darkorange", "blue"))+
  scale_fill_manual(values = c("darkorange", "blue"))+          
  labs(x = "",
       y = "Module eigengene")+  theme_classic()+ ggtitle("blue: 1259 genes")+
  theme(legend.position = "none", text =element_text(size =10))

ggsave("imgs/blue_male_module.png", t3, width = 2.5, height = 2.5)



t4 <- ME_df %>% filter(Module == "red") %>% 
  ggplot(aes(treatment, value, fill = treatment, color = treatment))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c( "darkorange", "blue"))+
  scale_fill_manual(values = c("darkorange", "blue"))+           
  labs(x = "",
       y = "Module eigengene")+  theme_classic()+ ggtitle("red: 713 genes")+
  theme(legend.position = "none", text =element_text(size =10))

ggsave("imgs/red_male_module.png", t4, width = 2.5, height = 2.5)




p1 <- ME_df %>% 
  ggplot(aes(con_ng_ml, value))+
  geom_point()+
  geom_smooth(method = "lm", se = F)+
  # scale_color_manual(values = c("blue", "darkorange"))+
  # scale_fill_manual(values = c("blue", "darkorange"))+         
  labs(x = "",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol = 5) + theme_classic()+
  theme(legend.position = "none", text =element_text(size =10))

p1


###########linear models 
head(ME_df)

ME_df$Module <- gsub("ME", "", ME_df$Module)


#blue
bl <- ME_df %>% filter(Module == "blue")

mb <- lm(value~treatment +con_ng_ml, data = bl)
summary(mb)
# treatmentNMX -0.17196    0.08851  -1.943   0.0663 .


#brown
br <- ME_df %>% filter(Module == "brown")

mbr <- lm(value~treatment +con_ng_ml, data = br)
summary(mbr)
# treatmentNMX  0.18425    0.08174   2.254    0.035 *
# #with T
# treatmentNMX  0.16831    0.08823   1.908   0.0709 .
# con_ng_ml     0.02642    0.04893   0.540   0.5951 


#turquoise
tur <- ME_df %>% filter(Module == "turquoise")

mt <- lm(value~treatment +con_ng_ml, data = tur)
summary(mt)
# treatmentNMX -0.19495    0.08834  -2.207   0.0392 *
#   con_ng_ml     0.01210    0.04899   0.247   0.8074 


#magenta
mag <- ME_df %>% filter(Module == "magenta")

mm <- lm(value~treatment+con_ng_ml, data = mag)
summary(mm)
# treatmentNMX  0.15634    0.08920   1.753   0.0950 .
# con_ng_ml    -0.08776    0.04946  -1.774   0.0912 .

#GREEN
g<- ME_df %>% filter(Module == "green")

mg <- lm(value~treatment+con_ng_ml, data = g)
summary(mg)
#nonsig 

#pink
p<- ME_df %>% filter(Module == "pink")

mp <- lm(value~treatment+con_ng_ml, data = p)
summary(mp)
# non sig

#red
r<- ME_df %>% filter(Module == "red")

mr <- lm(value~treatment+con_ng_ml, data = r)
summary(mr)
# treatmentNMX  0.20099    0.08157   2.464   0.0229 *
# con_ng_ml     0.04149    0.04523   0.917   0.3700  


#yellow
y<- ME_df %>% filter(Module == "yellow")

my <- lm(value~treatment+con_ng_ml, data = y)
summary(my)
#nonsig

#black
bla<- ME_df %>% filter(Module == "black")

mbla <- lm(value~treatment+con_ng_ml, data = bla)
summary(mbla)
#nonsig

#########module number. 


datExpr <- readRDS("males/clean_data/WGCNA/WGCNA_datExpr_males.RDS") 
net <- readRDS("males/clean_data/WGCNA/WGCNA_net_males_Power4.RDS")

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

moduleNumber = length(unique(moduleColors))

modNames = substring(names(MEs), 3)


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



#PART 5 =====================================================================================

# Module of interest via significance 
my_trait = "Status"
module = "brown"
module = "blue"
module = 'turquoise'
module = 'red'

module_list = c("brown",  "blue", "turquoise", "red")
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


write.csv(hubgenes_df, "males/clean_data/WGCNA/WCGNA_hubgene_list_TES_Power4.csv", row.names = F)


ADJ=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ, moduleColors)

Alldegrees1 %>%
  rownames_to_column('gene')-> kIN_df

xx <- cbind(gene = colnames(datExpr), module = moduleColors) %>%
  as.tibble()


kIN_df %>%
  left_join(xx) -> kIN_df

saveRDS(kIN_df,"males/clean_data/kIN_dataframe_TES_power4.RDS")


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
          "males/clean_data/WGCNA/WGCNA_all_gos_catogeryBP_TES_Power4.csv",
          row.names = F)


### hubgenes
library(tidyverse)
my_logFC_threshold = 0.2

# tes
tes<- readRDS("males/clean_data/limma_eFDR_Testis.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val) %>% distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold) %>%
  filter(.,P.Value <0.05)


t_up <- tes %>% filter(logFC >= 0.2) %>% arrange(-logFC) #65
t_down <- tes %>% filter(logFC <= -0.2)%>% arrange(logFC) #88


hub <- read_csv("males/clean_data/WGCNA/WCGNA_hubgene_list_TES_Power4.csv")
head(hub)

hbr <- hub %>% filter(moduleName == "brown") %>% arrange(-moduleMembership) #60 out of 840
head(hbr, 10)
hbr$gene[hbr$gene %in% t_down$genes] #none
hbr$gene[hbr$gene %in% t_up$genes]#10 

hb <- hub %>% filter(moduleName == "blue") %>% arrange(-moduleMembership) #125 out of 1259
head(hb, 10)
hb$gene[hb$gene %in% t_down$genes] #[1] "Eri3"
hb$gene[hb$gene %in% t_up$genes] #0
hb %>% filter(gene == "Eri3")

ht <- hub %>% filter(moduleName == "turquoise") %>% arrange(-moduleMembership) #85 out of 1301
head(ht, 10)

ht$gene[ht$gene %in% t_down$genes] #Hspb9
ht$gene[ht$gene %in% t_up$genes] #0
t_down %>% filter(genes %in% c("Eri3", "Hspb9"))

hs <- hub %>% filter(moduleName == "red") %>% arrange(-moduleMembership) #36 out of 713
head(hs, 10)
hs$gene[hs$gene %in% t_down$genes] #0
hs$gene[hs$gene %in% t_up$genes]#20


#look at overlap between brain regions and ovaries just female data. 
library(limma)
library(DESeq2)
library(edgeR)
library(AnnotationDbi)
library(tidyverse)
organism = 'org.Rn.eg.db'
library('org.Rn.eg.db')

my_logFC_threshold = 0.5
#ARC
arc_list <- readRDS("females/clean_data/limma_females_ARC_puberty.RDS")%>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~select(.,1:3,5))
ag <-arc_list$group
ap <-arc_list$puberty

ag$logFC <- ag$logFC*-1
ag_up <- ag %>% filter(logFC > 0.2)
ag_down <- ag %>% filter(logFC < 0.2)


ap_up <- ap %>% filter(logFC > 0.2)
ap_down <- ap %>% filter(logFC < 0.2)


#VMN
vmn_list <- readRDS("females/clean_data/limma_females_VMN_puberty.RDS")%>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~select(.,1:3,5))

vg <-vmn_list$group
vp <-vmn_list$puberty

vg$logFC <- vg$logFC*-1
vg_up <- vg %>% filter(logFC > 0.2)
vg_down <- vg %>% filter(logFC < 0.2)


vp_up <- vp %>% filter(logFC > 0.2)
vp_down <- vp %>% filter(logFC < 0.2)


#OVT
ovt_list <- readRDS("females/clean_data/limma_females_OVT_puberty.RDS")%>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~select(.,1:3,5))


og <-ovt_list$group
op <-ovt_list$puberty

og$logFC <- og$logFC*-1
og_up <- og %>% filter(logFC > 0.2)
og_down <- og %>% filter(logFC < 0.2)


op_up <- op %>% filter(logFC > 0.2)
op_down <- op %>% filter(logFC < 0.2)


#get genes for boxplots: 
vmn_gp_down <- vg_down[vg_down$genes %in% vp_down$genes,] %>% select(genes)


#get expression data 
ex <- readRDS("females/clean_data/limma_vdl_females_VMN_puberty.RDS")
head(ex)

x <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('gene') 

ide<- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids)

# sample information 
id <- read_csv("raw_data/traitData.csv")
colnames(id)
idx <- id %>% select(1:4,15)
v_data <- idx %>% filter(Sex =="F")
v_data$SampleName <- paste("V", v_data$SampleName, sep = "_")

v_datax <- v_data %>% filter(SampleName %in% ide$ids)

x %>% 
  filter(gene %in% vmn_gp_down$genes) -> xex
colnames(xex)

xex2 <- xex %>% pivot_longer(cols = 2:18, names_to = "SampleName")

p <- xex2 %>% full_join(v_datax)

p$group <- factor(p$group, levels = c("Blue", "Orange"))

source("functions/geom_boxjitter.R")
p1 <- ggplot(p, aes(Group,value, color = Group, fill = Group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA, alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE, position = position_dodge(0.85)) +
  scale_color_manual(values = c("blue", "darkorange"))+
  scale_fill_manual(values = c("blue", "darkorange"))+
  facet_wrap(~gene, ncol = 6)+
scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none", text =element_text(size = 15))
p1



p1 <- ggplot(p, aes(Puberty,value))+
  geom_point(color = "blue", fill = "blue", size = 1) +
  stat_smooth(method = "lm", se = F)+
  facet_wrap(~gene, ncol = 6)+
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none", text =element_text(size = 15))
p1

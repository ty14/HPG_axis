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
arc_list <- readRDS("males/clean_data/limma_males_ARC_puberty.RDS")%>% 
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
vmn_list <- readRDS("males/clean_data/limma_males_VMN_puberty.RDS")%>% 
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


#TES
Tes_list <- readRDS("males/clean_data/limma_males_Testis_puberty.RDS")%>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~select(.,1:3,5))


tg <-Tes_list$group
tp <-Tes_list$puberty

tg$logFC <- tg$logFC*-1
tg_up <- tg %>% filter(logFC > 0.2)
tg_down <- tg %>% filter(logFC < 0.2)


tp_up <- tp %>% filter(logFC > 0.2)
tp_down <- tp %>% filter(logFC < 0.2)

# ARC and Tes group overlap # 1 overlap down in Tes up in ARC  

# upsetter Plot 
library(UpSetR)
library(workflowr)
library(ComplexUpset)



##########now look between group and puberty

#TES - genes opposite directions
listInput <- list(tg_up= tg_up$genes, tp_up=tp_up$genes,
                  tg_down = tg_down$genes, tp_down = tp_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)


tg_down[tg_down$genes %in% tp_up$genes,]


#VMN
listInput <- list(vg_up= vg_up$genes, vp_up=vp_up$genes,
                  vg_down = vg_down$genes, vp_down = vp_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)


# gene opposite direction. 
vg_up[vg_up$genes %in% vp_down$genes,]


#ARC
listInput <- list(ag_up= ag_up$genes, ap_up=ap_up$genes,
                  ag_down = ag_down$genes, ap_down = ap_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

ap_up[ap_up$genes %in% ag_down$genes,]


# VMN and ARC group overlap #2 genes same direction.  

listInput <- list(ag_up= ag_up$genes, vg_up=vg_up$genes,
                  ag_down = ag_down$genes, vg_down = vg_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)
ag_down[ag_down$genes %in% vg_down$genes,]


# VMN and ARC puberty overlap # 2 genes

listInput <- list(ap_up= ap_up$genes, vp_up=vp_up$genes,
                  ap_down = ap_down$genes, vp_down = vp_down$genes)

UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

ap_down[ap_down$genes %in% vp_down$genes,]


#VMN and Testis # 2 genes 
listInput <- list(tg_up= tg_up$genes, ag_up=ag_up$genes,
                  tg_down = tg_down$genes, ag_down = ag_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

ag_down[ag_down$genes %in% tg_down$genes,]

# VMN and Tes puberty overlap # 1 genes

listInput <- list(tp_up= tp_up$genes, vp_up=vp_up$genes,
                  tp_down = tp_down$genes, vp_down = vp_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

vp_down[vp_down$genes %in% tp_down$genes,]

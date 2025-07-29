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



#ARC
arc_list <- readRDS("females/clean_data/limma_eFDR_ARC.RDS")%>% 
  distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold) %>%
  filter(.,P.Value <0.05) %>% 
  select(.,1:3,5)
ag1 <-arc_list


ag_up1 <- ag1 %>% filter(logFC > 0.2)
ag_down1 <- ag1 %>% filter(logFC < 0.2)



#VMN
vmn_list <- readRDS("females/clean_data/limma_eFDR_VMN.RDS")%>% 
  distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold) %>%
  filter(.,P.Value <0.05) %>% 
  select(.,1:3,5)

vg1 <-vmn_list


vg_up1 <- vg1 %>% filter(logFC > 0.2)
vg_down1 <- vg1 %>% filter(logFC < 0.2)


#OVT
ovt_list <- readRDS("females/clean_data/limma_eFDR_OVT.RDS")%>% 
  distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold) %>%
  filter(.,P.Value <0.05) %>% 
  select(.,1:3,5)

og1 <-ovt_list


og_up1 <- og1 %>% filter(logFC > 0.2)
og_down1 <- og1 %>% filter(logFC < 0.2)



## LOOK AT OVERLAP IN TWO ANALYSIS AT EACH REGION. 

# OVT
# upsetter Plot 
library(UpSetR)
library(workflowr)
library(ComplexUpset)

listInput <- list(og_up= og_up$genes, og_up1=og_up1$genes,
                  og_down = og_down$genes, og_down1 = og_down1$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)



# upsetter Plot for VMN
listInput <- list(vg_up= vg_up$genes, vg_up1=vg_up1$genes,
                  vg_down = vg_down$genes, vg_down1 = vg_down1$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)


# upsetter Plot for ARC
listInput <- list(ag_up= ag_up$genes, ag_up1=ag_up1$genes,
                  ag_down = ag_down$genes, ag_down1 = ag_down1$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)



##########now look between group and puberty

#OVT
listInput <- list(og_up= og_up$genes, op_up=op_up$genes,
                  og_down = og_down$genes, op_down = op_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

og_down[og_down$genes %in% op_up$genes,]


#VMN
listInput <- list(vg_up= vg_up$genes, vp_up=vp_up$genes,
                  vg_down = vg_down$genes, vp_down = vp_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

vmn_gp_down <- vg_down[vg_down$genes %in% vp_down$genes,]

vg_down[vg_down$genes %in% vp_up$genes,]
vg_up[vg_up$genes %in% vp_down$genes,]


#ARC
listInput <- list(ag_up= ag_up$genes, ap_up=ap_up$genes,
                  ag_down = ag_down$genes, ap_down = ap_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

ap_up[ap_up$genes %in% ag_down$genes,]


# VMN and ARC group overlap # none 

listInput <- list(ag_up= ag_up$genes, vg_up=vg_up$genes,
                  ag_down = ag_down$genes, vg_down = vg_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)


# VMN and ARC puberty overlap # 7 genes

liSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)


ap_up[ap_up$genes %in% vp_up$genes,]

ap_down[ap_down$genes %in% vp_down$genes,]


# VMN and OVT group overlap # 5 overlap down in VMN up in ovary  

listInput <- list(og_up= og_up$genes, vg_up=vg_up$genes,
                  og_down = og_down$genes, vg_down = vg_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

og_up[og_up$genes %in% vg_down$genes,]

vg_down %>% filter(genes %in% c("Cdhr4", "Ccdc78", "Cfap99", "Lrrc34", "Ccdc180"))

# VMN and OVT puberty overlap # none

listInput <- list(op_up= op_up$genes, vp_up=vp_up$genes,
                  op_down = op_down$genes, vp_down = vp_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)
stInput <- list(ap_up= ap_up$genes, vp_up=vp_up$genes,
                  ap_down = ap_down$genes, vp_down = vp_down$genes)


Up


# ARC and OVT group overlap # 1 overlap down in OVT up in ARC  

listInput <- list(og_up= og_up$genes, ag_up=ag_up$genes,
                  og_down = og_down$genes, ag_down = ag_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

ag_up[ag_up$genes %in% og_down$genes,]

# VMN and OVT puberty overlap #  genes

listInput <- list(op_up= op_up$genes, ap_up=ap_up$genes,
                  op_down = op_down$genes, ap_down = ap_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

ap_up[ap_up$genes %in% op_up$genes,]
ap_down[ap_down$genes %in% op_down$genes,]

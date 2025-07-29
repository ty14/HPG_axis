
# go analysis
library(clusterProfiler)
library(AnnotationDbi)
library(org.Rn.eg.db)
library(tidyverse)

#Descenders 
my_logFC_threshold = 0.2

y1a <- readRDS("females/clean_data/limma_eFDR_ARC.RDS") %>%
  distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold) 

y2a <- readRDS("females/clean_data/limma_eFDR_VMN.RDS") %>% 
  distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold) 

y3a <- readRDS("females/clean_data/limma_eFDR_OVT.RDS")%>% 
  distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold) 

y4a <- readRDS("females/clean_data/limma_eFDR_PIT.RDS")%>% 
  distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold)

y1a <- y1a %>% 
  mutate(entrez = mapIds(org.Rn.eg.db, keys = y1a$genes,
                         column = "ENTREZID", keytype = "SYMBOL")) %>% 
  filter(!is.na(entrez))

y2a <- y2a %>% 
  mutate(entrez = mapIds(org.Rn.eg.db, keys = y2a$genes,
                         column = "ENTREZID", keytype = "SYMBOL")) %>% 
  filter(!is.na(entrez))

y3a <- y3a %>% 
  mutate(entrez = mapIds(org.Rn.eg.db, keys = y3a$genes,
                         column = "ENTREZID", keytype = "SYMBOL")) %>% 
  filter(!is.na(entrez))

y4a <- y4a %>% 
  mutate(entrez = mapIds(org.Rn.eg.db, keys = y4a$genes,
                         column = "ENTREZID", keytype = "SYMBOL")) %>% 
  filter(!is.na(entrez))




source("functions/gettop10GO_Rat.R")

gettop10GO(y1a, my_showCategory) %>% 
  mutate(tissue = "ARC") -> top10go1

gettop10GO(y2a, my_showCategory) %>% 
  mutate(tissue = "VMN") -> top10go2

gettop10GO(y3a, my_showCategory) %>% 
  mutate(tissue = "OVT") -> top10go3

gettop10GO(y4a, my_showCategory) %>% 
  mutate(tissue = "PIT") -> top10go4


top10go_all <- top10go1 %>% rbind(top10go2, top10go3, top10go4)

write.csv(top10go_all, "females/clean_data/GOTerms_females_limma_allgenes_LogFC0.2.csv", row.names = F)



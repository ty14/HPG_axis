library(tidyverse)
# sample information 

df <- read_csv("females/clean_data/RIA_estrodial.csv")
head(df)

colnames(df)<- c("SampleName", "treatment", "mean", "sd", "cv", "con_pg_ml")
head(df)

id <- read_csv("raw_data/traitData.csv")
colnames(id)
idx <- id %>% dplyr::select(1:4,15)


dfx <- df %>% full_join(idx)


fe_data <- dfx %>% dplyr::select(1,2,6,7:10)

e2_genes <- c("Kiss1", "Tac3", "Pdyn", "Gal", "Lepr", "Ghrh", "Ghr", "Esr1", "Esr2", "Pomc", "Ar", "Nr5a2", 
  "Ghrh", "Ghr", "Esr1", "Esr2","Cyp19a1", "Cyp17a1", "Cyp11a1", "Hsd17b1", "Star", "Fshb", "Fshr",
  "Lhb", "Lhcgr","Amh", "Inha",  "Ldlr", "Scarb1", "Hmgcr")

# look for hubgenes in modules that correlate with E2

fe_hub <- read.csv("females/clean_data/WGCNA/WCGNA_hubgene_list_OVT_Power4.csv")
head(fe_hub)


#modules to look at 
# black
# brown
# cyan 
# darkgreen
# lightgreen 


bh <- fe_hub %>% filter(moduleName == "black", moduleMembership > 0.90)
#42 genes

brh <- fe_hub %>% filter(moduleName == "brown", moduleMembership > 0.90)
#188 genes

dgh <- fe_hub %>% filter(moduleName == "darkgreen", moduleMembership > 0.80)
#5 genes

#all unique 

all_e2_genes <- c(e2_genes, bh$gene, brh$gene, dgh$gene)



#first ARC
#ARC
arc_list <- readRDS("females/clean_data/limma_vdl_ARC.RDS")


x <- arc_list$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('gene') 

ide<- arc_list$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids)

x %>% 
  filter(gene %in% e2_genes) -> xex 
colnames(xex)
dim(xex) # 49 GENES 
xex2 <- xex %>% pivot_longer(cols = 2:12, names_to = "SampleName")
xex2$SampleName<- gsub("A_", "",xex2$SampleName )

p <- xex2 %>% full_join(fe_data)
p <- na.omit(p)

ARC_e2 <- ggplot(p, aes(con_pg_ml,value))+
  geom_point(color = "blue", fill = "blue", size = 1) +
  stat_smooth(method = "lm", se = F)+
  facet_wrap(~gene, ncol = 5) +
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("Estradiol pg/ml")+
  theme_bw()+
  ggtitle("ARC")+
  theme(legend.position = "none", text =element_text(size = 15))

ggsave("imgs/ARC_e2_use.png", ARC_e2, width = 8, height = 5.5)


tac <- p %>% filter(gene == "Tac3")
cor.test(tac$con_pg_ml, tac$value)


k <- p %>% filter(gene == "Kiss1")
cor.test(k$con_pg_ml, k$value)

lh <- p %>% filter(gene == "Lhb")
cor.test(lh$con_pg_ml, lh$value)


#### VMN
vmn_list <- readRDS("females/clean_data/limma_vdl_VMN.RDS")


x <- vmn_list$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('gene') 

ide<- vmn_list$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids)

x %>% 
  filter(gene %in%e2_genes) -> xex
colnames(xex)
dim(xex)
# 51 genes

xex2 <- xex %>% pivot_longer(cols = 2:18, names_to = "SampleName")
xex2$SampleName<- gsub("V_", "",xex2$SampleName )

p <- xex2 %>% full_join(fe_data)
p <- na.omit(p)

vmn_e2 <- ggplot(p, aes(con_pg_ml,value))+
  geom_point(color = "blue", fill = "blue", size = 1) +
  stat_smooth(method = "lm", se = F)+
  facet_wrap(~gene, ncol = 5)+
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("Estradiol pg/ml")+
  theme_bw()+
  ggtitle("VMN")+
  theme(legend.position = "none", text =element_text(size = 15))


ggsave("imgs/VMN_e2_use.png", vmn_e2, width = 8, height = 5.5)

##
gr <- p %>% filter(gene == "Ghrh")
cor.test(gr$con_pg_ml, gr$value)


lh <- p %>% filter(gene == "Lhb")
cor.test(lh$con_pg_ml, lh$value)

po <- p %>% filter(gene == "Pomc")
cor.test(po$con_pg_ml, po$value)

lep <- p %>% filter(gene == "Lepr")
cor.test(lep$con_pg_ml, lep$value)



#Pit
pit_list <- readRDS("females/clean_data/limma_vdl_PIT.RDS")


x <- pit_list$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('gene') 

ide<- pit_list$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids)

x %>% 
  filter(gene %in% e2_genes) -> xex
colnames(xex)

xex2 <- xex %>% pivot_longer(cols = 2:24, names_to = "SampleName")
xex2$SampleName<- gsub("V_", "",xex2$SampleName )
dim(xex)
p <- xex2 %>% full_join(fe_data)
p <- na.omit(p)

PIT_e2 <- ggplot(p, aes(con_pg_ml,value))+
  geom_point(color = "blue", fill = "blue", size = 1) +
  stat_smooth(method = "lm", se = F)+
  facet_wrap(~gene, ncol = 5)+
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("Estradiol pg/ml")+
  theme_bw()+
  ggtitle("Pituitary")+
  theme(legend.position = "none", text =element_text(size = 15))

ggsave("imgs/PIT_e2_use.png", PIT_e2, width = 8, height = 3.75)


gr <- p %>% filter(gene == "Ghr")
cor.test(gr$con_pg_ml, gr$value)


gal <- p %>% filter(gene == "Gal")
cor.test(gal$con_pg_ml, gal$value)

fsh <- p %>% filter(gene == "Fshb")
cor.test(fsh$con_pg_ml, fsh$value)



##OVary 

#OV
OV_list <- readRDS("females/clean_data/limma_vdl_OVT.RDS")


x <- OV_list$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('gene') 

ide<- OV_list$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids)

x %>% 
  filter(gene %in% e2_genes) -> xex
colnames(xex)

xex2 <- xex %>% pivot_longer(cols = 2:24, names_to = "SampleName")
xex2$SampleName<- gsub("V_", "",xex2$SampleName )
dim(xex)
p <- xex2 %>% full_join(fe_data)
p <- na.omit(p)

OV_e2 <- ggplot(p, aes(con_pg_ml,value))+
  geom_point(color = "blue", fill = "blue", size = 1) +
  stat_smooth(method = "lm", se = F)+
  facet_wrap(~gene, ncol = 5)+
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("Estradiol pg/ml")+
  theme_bw()+
  ggtitle("Ovary")+
  theme(legend.position = "none", text =element_text(size = 15))

ggsave("imgs/OVT_e2_use.png", OV_e2, width = 8, height = 5.5)



library(tidyverse)

  df <- read_csv("females/clean_data/RIA_estrodial.csv")
  head(df)
  
  colnames(df)<- c("SampleName", "treatment", "mean", "sd", "cv", "con_pg_ml")
  head(df)

df2<- df %>% mutate(sd_con = sd(con_pg_ml, na.rm = T))

source("functions/geom_boxjitter.R")

p <- df2 %>%  ggplot(aes(treatment, con_pg_ml, fill = treatment, color = treatment))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("darkorange", "blue"))+
  scale_fill_manual(values = c("darkorange", "blue"))+         
  labs(x = "",
       y = "Estradiol pg/ml")+
  theme_classic()+
  theme(legend.position = "none", text =element_text(size =15))
  
p
ggsave("imgs/estradiol_boxplot.png", p, width = 3, height = 3,  dpi = 300)

library(lmerTest)

hist(df$con_pg_ml)

m1 <- lmer(con_pg_ml ~ treatment+ (1|Cohort), data = fe_data)
summary(m1)


df3 <- df2 %>% pivot_wider(names_from = "treatment", values_from = con_pg_ml)
t.test(df3$Control, df3$NMX, na.rm =T, var.equal = T)



# sample information 
id <- read_csv("raw_data/traitData.csv")
colnames(id)
idx <- id %>% dplyr::select(1:4,15)


dfx <- df %>% full_join(idx)


fe_data <- dfx %>% dplyr::select(1,2,6,7:10)

write.csv(fe_data, "females/clean_data/female_trait_data.csv", row.names = F)

#estrdiol range and mean 

range(dfx$con_pg_ml, na.rm = T)
mean(dfx$cv, na.rm = T)

median(df3$Control, na.rm = T)
median(df3$NMX, na.rm = T)

quantile(df3$Control, 0.25,na.rm = T)
quantile(df3$Control, 0.75,na.rm = T)

quantile(df3$NMX, 0.25,na.rm = T)
quantile(df3$NMX, 0.75,na.rm = T)

mean(df3$NMX, na.rm = T)
sd(df3$NMX, na.rm = T)
length(df3$NMX)
sd(df3$NMX, na.rm = T)/12

mean(df3$Control, na.rm = T)
sd(df3$Control, na.rm = T)
unique(df3$Control)
sd(df3$Control, na.rm = T)/11


##genes that correlate with E2?

library(tidyverse)
my_logFC_threshold = 0.2
  
# females vmn
fe_vmn <- readRDS("females/clean_data/limma_eFDR_VMN.RDS") %>% 
  select(genes, logFC, P.Value,  adj.P.Val)  %>% distinct(.) %>% 
filter(.,abs(logFC) >= my_logFC_threshold)
# females arc
fe_arc <- readRDS("females/clean_data/limma_eFDR_ARC.RDS") %>% 
  select(genes, logFC, P.Value,  adj.P.Val) %>% distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold)

# females pit
fe_pit <- readRDS("females/clean_data/limma_eFDR_PIT.RDS") %>% 
  select(genes, logFC, P.Value,  adj.P.Val)%>% distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold)


fe_ag <- fe_pit %>% filter(genes %in% c("Lepr1", "Tac3", "Esr1", "Gal", "Lepr", "Pdyn", "Npy", "Oxt", "Ntsr1", "Ntsr2", "Pgr", "Nr5a1", "Bdnf", "Lh", "Fsh", "Cga","Prl", "Pomc", "Gh", "Pou1f1", "Gnrhr", "Avp"))
fe_data



# looking at correlation of genes with E2

my_logFC_threshold = 0.2
#ARC
arc_list <- readRDS("females/clean_data/limma_vdl_ARC.RDS")


x <- arc_list$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('gene') 

ide<- arc_list$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids)

x %>% 
  filter(gene %in% c("Kiss1", "Tac3", "Pdyn", "Gal", "Lepr", "Ghrh", "Ghr", "Esr1", "Esr2", "Pomc", "Ar", "Nr5a2")) -> xex
colnames(xex)

xex2 <- xex %>% pivot_longer(cols = 2:13, names_to = "SampleName")
xex2$SampleName<- gsub("A_", "",xex2$SampleName )

p <- xex2 %>% full_join(fe_data)
p <- na.omit(p)

ggplot(p, aes(con_pg_ml,value))+
  geom_point(color = "blue", fill = "blue", size = 1) +
  stat_smooth(method = "lm", se = F)+
  facet_wrap(~gene, ncol = 4)+
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("Estradiol pg/ml")+
  theme_bw()+
  theme(legend.position = "none", text =element_text(size = 15))


# stats
kiss1 <- p %>% filter(gene == "Kiss1")

library(lmerTest)
km <- lm(value~con_pg_ml, data = kiss1)
summary(km)

cor.test(kiss1$value, kiss1$con_pg_ml)

lepr <- p %>% filter(gene == "Lepr")

library(lmerTest)
lm <- lm(value~con_pg_ml, data = lepr)
summary(lm)

cor.test(lepr$value, lepr$con_pg_ml)

#vmn
#ARC
vmn_list <- readRDS("females/clean_data/limma_vdl_VMN.RDS")


x <- vmn_list$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('gene') 

ide<- vmn_list$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids)

x %>% 
  filter(gene %in% c("Kiss1", "Tac3", "Pdyn", "Gal", "Lepr", "Ghrh", "Ghr", "Esr1", "Esr2", "Pomc", "Ar", "Nr5a2")) -> xex
colnames(xex)

xex2 <- xex %>% pivot_longer(cols = 2:18, names_to = "SampleName")
xex2$SampleName<- gsub("V_", "",xex2$SampleName )

p <- xex2 %>% full_join(fe_data)
p <- na.omit(p)

ggplot(p, aes(con_pg_ml,value))+
  geom_point(color = "blue", fill = "blue", size = 1) +
  stat_smooth(method = "lm", se = F)+
  facet_wrap(~gene, ncol = 4)+
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("Estradiol pg/ml")+
  theme_bw()+
  theme(legend.position = "none", text =element_text(size = 15))


# stats
pomc <- p %>% filter(gene == "Lepr")

library(lmerTest)
km <- lm(value~con_pg_ml, data = pomc)
summary(km)

cor.test(pomc$value, pomc$con_pg_ml)

lepr <- p %>% filter(gene == "Lepr")


#Pit
pit_list <- readRDS("females/clean_data/limma_vdl_PIT.RDS")


x <- pit_list$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('gene') 

ide<- pit_list$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids)

x %>% 
  filter(gene %in% c("Kiss1", "Tac3", "Pdyn", "Gal", "Lepr", "Ghrh", "Ghr", "Esr1", "Esr2", "Pomc", "Ar", "Nr5a2")) -> xex
colnames(xex)

xex2 <- xex %>% pivot_longer(cols = 2:24, names_to = "SampleName")
xex2$SampleName<- gsub("V_", "",xex2$SampleName )

p <- xex2 %>% full_join(fe_data)
p <- na.omit(p)

ggplot(p, aes(con_pg_ml,value))+
  geom_point(color = "blue", fill = "blue", size = 1) +
  stat_smooth(method = "lm", se = F)+
  facet_wrap(~gene, ncol = 6)+
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("Estradiol pg/ml")+
  theme_bw()+
  theme(legend.position = "none", text =element_text(size = 15))


# stats
pomc <- p %>% filter(gene == "Gal")

library(lmerTest)
km <- lm(value~con_pg_ml, data = pomc)
summary(km)

cor.test(pomc$value, pomc$con_pg_ml)

lepr <- p %>% filter(gene == "Lepr")

library(lmerTest)
lm <- lm(value~con_pg_ml, data = lepr)
summary(lm)

cor.test(lepr$value, lepr$con_pg_ml)


#OVT
ovt_list <- readRDS("females/clean_data/limma_vdl_OVT.RDS")


x <- ovt_list$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('gene') 

ide<- ovt_list$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids)

x %>% 
  filter(gene %in% c("Ghrh", "Ghr", "Esr1", "Esr2","Cyp19a1", "Cyp17a1", "Cyp11a1", "Hsd17b1", "Star", "Fshb", "Fshr", "Lh", "Lhcgr","Amh", "Inha",  "Ldlr", "Scarb1", "Hmgcr")) -> xex
colnames(xex)

xex2 <- xex %>% pivot_longer(cols = 2:24, names_to = "SampleName")
xex2$SampleName<- gsub("V_", "",xex2$SampleName )

p <- xex2 %>% full_join(fe_data)
p <- na.omit(p)

ggplot(p, aes(con_pg_ml,value))+
  geom_point(color = "blue", fill = "blue", size = 1) +
  stat_smooth(method = "lm", se = F)+
  facet_wrap(~gene, ncol = 6)+
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("Estradiol pg/ml")+
  theme_bw()+
  theme(legend.position = "none", text =element_text(size = 15))




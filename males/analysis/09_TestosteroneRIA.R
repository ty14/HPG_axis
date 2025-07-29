library(tidyverse)

source("functions/geom_boxjitter.R")

library(readxl)
df<- read_excel("raw_data/250206_TestosteroneRIA 1.xlsx")
head(df)
df <- df[,c(1:5,8)]
colnames(df)<- c("SampleName", "treatment", "mean", "sd", "cv", "con_ng_ml")
head(df)

df2<- df %>% mutate(sd_con = sd(con_ng_ml, na.rm = T))

df$SampleName


nmx <- df2[c(74:97),]
nmx$SampleName
mean(nmx$cv, na.rm = T) #8.8

nmx2 <- nmx %>% filter(SampleName != "NMX23B")
nmx2$SampleName

mean(nmx2$cv, na.rm = T) #8.0


# sample information 
id <- read_csv("raw_data/traitData.csv")
colnames(id)
idx <- id %>% dplyr::select(1:4,15)

#remove NMX08A  for being an outlier 
dfx <- nmx2 %>% full_join(idx) %>% filter(SampleName != "NMX08A")%>% 
  filter(SampleName != "NMX23B" )%>% filter(Sex != "F")


m_data <- dfx %>% dplyr::select(1,2,6,7:10)
table(m_data$treatment)

write.csv(m_data, "males/clean_data/male_trait_data.csv", row.names = F)
m_data

p <- m_data %>%  ggplot(aes(treatment, con_ng_ml, fill = treatment, color = treatment))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("darkorange", "blue"))+
  scale_fill_manual(values = c("darkorange", "blue"))+         
  labs(x = "",
       y = "Testosterone ng/ml")+
  theme_classic()+
  theme(legend.position = "none", text =element_text(size =15))

p
ggsave("imgs/test_boxplot.png", p, width = 3, height = 3, dpi = 300)



library(lmerTest)
hist(m_data$con_ng_ml)

m1 <- lmer(con_ng_ml ~ treatment+ (1|Cohort), data = m_data)
summary(m1)

c <- m_data %>% filter(treatment == "Control")

median(c$con_ng_ml)
quantile(c$con_ng_ml, 0.25)
quantile(c$con_ng_ml, 0.75)

n <- m_data %>% filter(treatment == "NMX")

median(n$con_ng_ml)
quantile(n$con_ng_ml, 0.25)
quantile(n$con_ng_ml, 0.75)



library(tidyverse)
# sample information 

df <- read_csv("males/clean_data/male_trait_data.csv")


T_genes <- c("Kiss1", "Tac3", "Pdyn", "Gal", "Lepr", "Ghrh", "Ghr", "Esr1", "Esr2", "Pomc", "Ar", "Nr5a2", 
              "Ghrh", "Ghr", "Esr1", "Esr2","Cyp19a1", "Cyp17a1", "Cyp11a1", "Hsd17b1", "Star", "Fshb", "Fshr",
              "Lhb", "Lhcgr","Amh", "Inha",  "Ldlr", "Scarb1", "Hmgcr")



#ARC
arc_list <- readRDS("males/clean_data/limma_vdl_ARC.RDS")


x <- arc_list$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('gene') 

ide<- arc_list$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids)

x %>% 
  filter(gene %in% T_genes) -> xex 
colnames(xex)
dim(xex) # 49 GENES 
xex2 <- xex %>% pivot_longer(cols = 2:14, names_to = "SampleName")
xex2$SampleName<- gsub("A_", "",xex2$SampleName )

p <- xex2 %>% full_join(df)
p <- na.omit(p)

ARC_t <- ggplot(p, aes(con_ng_ml,value))+
  geom_point(color = "blue", fill = "blue", size = 1) +
  stat_smooth(method = "lm", se = F)+
  facet_wrap(~gene, ncol = 5) +
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("Testosterone ng/ml")+
  theme_bw()+
  ggtitle("ARC")+
  theme(legend.position = "none", text =element_text(size = 15))

ARC_t

ggsave("imgs/ARC_tes_use.png", ARC_t, width = 8, height = 5.5)



ARC_t2 <- ggplot(p, aes(con_ng_ml,value, color = treatment))+
  geom_point( size = 1) +
  stat_smooth(method = "lm", se = F)+
  facet_wrap(~gene, ncol = 5) +
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("Testosterone ng/ml")+
  theme_bw()+
  ggtitle("ARC")+
  theme(legend.position = "none", text =element_text(size = 15))

ARC_t2

gh <- p %>% filter(gene == "Ghrh")
cor.test(gh$con_ng_ml, gh$value)
# #sig
# data:  gh$con_ng_ml and gh$value
# t = -5.4317, df = 8, p-value = 0.0006221
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.9731342 -0.5827837
# sample estimates:
#   cor 
# -0.8869545 
k <- p %>% filter(gene == "Lepr")
cor.test(k$con_ng_ml, k$value)

kc <- p %>% filter(gene == "Kiss1", treatment == "Control")
cor.test(kc$con_ng_ml, kc$value)



#vmn
vmn_list <- readRDS("males/clean_data/limma_vdl_VMN.RDS")


x <- vmn_list$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('gene') 

ide<- vmn_list$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids)

x %>% 
  filter(gene %in% T_genes) -> xex 
colnames(xex)
dim(xex) # 15 GENES 
xex2 <- xex %>% pivot_longer(cols = 2:20, names_to = "SampleName")
xex2$SampleName<- gsub("V_", "",xex2$SampleName )

p <- xex2 %>% full_join(df)
p <- na.omit(p)

VMN_t <- ggplot(p, aes(con_ng_ml,value))+
  geom_point(color = "blue", fill = "blue", size = 1) +
  stat_smooth(method = "lm", se = F)+
  facet_wrap(~gene, ncol = 5) +
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("Testosterone ng/ml")+
  theme_bw()+
  ggtitle("VMN")+
  theme(legend.position = "none", text =element_text(size = 15))

VMN_t

ggsave("imgs/VMN_tes_use.png", VMN_t, width = 8, height = 5.5)



VMN_t2 <- ggplot(p, aes(con_ng_ml,value, color = treatment))+
  geom_point( size = 1) +
  stat_smooth(method = "lm", se = F)+
  facet_wrap(~gene, ncol = 5) +
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("Testosterone ng/ml")+
  theme_bw()+
  ggtitle("VMN")+
  theme(legend.position = "none", text =element_text(size = 15))

VMN_t2


esr <- p %>% filter(gene == "Esr1")
cor.test(esr$con_ng_ml, esr$value)

k <- p %>% filter(gene == "Pdyn")
cor.test(k$con_ng_ml, k$value)


#pit
pit_list <- readRDS("males/clean_data/limma_vdl_PIT.RDS")


x <- pit_list$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('gene') 

ide<- pit_list$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids)

x %>% 
  filter(gene %in% T_genes) -> xex 
colnames(xex)
dim(xex) # 15 GENES 
xex2 <- xex %>% pivot_longer(cols = 2:25, names_to = "SampleName")
xex2$SampleName<- gsub("P_", "",xex2$SampleName )

p <- xex2 %>% full_join(df)
p <- na.omit(p)

pit_t <- ggplot(p, aes(con_ng_ml,value))+
  geom_point(color = "blue", fill = "blue", size = 1) +
  stat_smooth(method = "lm", se = F)+
  facet_wrap(~gene, ncol = 5) +
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("Testosterone ng/ml")+
  theme_bw()+
  ggtitle("Pituitary")+
  theme(legend.position = "none", text =element_text(size = 15))

pit_t

ggsave("imgs/pit_tes_use.png", pit_t, width = 8, height = 3.5)



pit_t2 <- ggplot(p, aes(con_ng_ml,value, color = treatment))+
  geom_point( size = 1) +
  stat_smooth(method = "lm", se = F)+
  facet_wrap(~gene, ncol = 5) +
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("Testosterone ng/ml")+
  theme_bw()+
  ggtitle("pit")+
  theme(legend.position = "none", text =element_text(size = 15))

pit_t2


esr <- p %>% filter(gene == "Ar", treatment == "Control")
cor.test(esr$con_ng_ml, esr$value)

k <- p %>% filter(gene == "Pomc")
cor.test(k$con_ng_ml, k$value)



#pit
tes_list <- readRDS("males/clean_data/limma_vdl_TES.RDS")


x <- tes_list$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('gene') 

ide<- tes_list$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids)

x %>% 
  filter(gene %in% T_genes) -> xex 
colnames(xex)
dim(xex) # 15 GENES 
xex2 <- xex %>% pivot_longer(cols = 2:24, names_to = "SampleName")
xex2$SampleName<- gsub("V_", "",xex2$SampleName )

p <- xex2 %>% full_join(df)
p <- na.omit(p)

tes_t <- ggplot(p, aes(con_ng_ml,value))+
  geom_point(color = "blue", fill = "blue", size = 1) +
  stat_smooth(method = "lm", se = F)+
  facet_wrap(~gene, ncol = 4) +
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("Testosterone ng/ml")+
  theme_bw()+
  ggtitle("Testis")+
  theme(legend.position = "none", text =element_text(size = 15))

tes_t

ggsave("imgs/Testis_tes_use.png", tes_t, width = 6, height = 3.5)



tes_t2 <- ggplot(p, aes(con_ng_ml,value, color = treatment))+
  geom_point( size = 1) +
  stat_smooth(method = "lm", se = F)+
  facet_wrap(~gene, ncol = 5) +
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("Testosterone ng/ml")+
  theme_bw()+
  ggtitle("TES")+
  theme(legend.position = "none", text =element_text(size = 15))

tes_t2


esr <- p %>% filter(gene == "Lhcgr")
cor.test(esr$con_ng_ml, esr$value)
#barely. 
k <- p %>% filter(gene == "Cyp11a1")
cor.test(k$con_ng_ml, k$value)


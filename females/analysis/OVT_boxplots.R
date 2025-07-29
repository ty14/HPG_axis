# genes in the OVT
# females 

ex <- readRDS("females/clean_data/limma_vdl_OVT.RDS")
head(ex)

x <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('gene') 


# sample information 
id <- read_csv("raw_data/traitData.csv")
colnames(id)
idx <- id %>% select(1:4,15)

#coldata arc
a_data <- idx
# a_data$SampleName <- paste("", a_data$SampleName, sep = "_")
a_data$ids <- a_data$SampleName
a_datax <- a_data %>% column_to_rownames(., var = "SampleName")
a_datax$bi_sex <- ifelse(a_datax$Sex == "M", 0, 1)
a_datax$bi_group <- ifelse(a_datax$Group == "Blue", 0, 1)


a_datax$Group



top_genes <- c("Cfap65", "Dnah12", "Dnah11", "Dnai7", "Celsr2", "Erich3")
down_genes <- c("RT1-DMa","RT1-DMb", "RT1-Db1", "Ccl7", "Coro1a", "Mpeg1")
x %>% 
  filter(gene %in% top_genes) -> xex
colnames(xex)

xex2 <- xex %>% pivot_longer(cols = 2:23, names_to = "ids")

p <- xex2 %>% full_join(a_data)

p$Group <- ifelse(p$Group == "Orange", "Control", "NMX")
p$group <- factor(p$Group, levels = c("Control", "NMX"))
p <- p %>% na.omit(.)

source('functions/geom_boxjitter.R')
library(viridis)

p1 <- ggplot(p, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA, alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE, position = position_dodge(0.85)) +
  scale_color_manual(values = c( "darkorange", "blue"))+
  scale_fill_manual(values = c( "darkorange", "blue"))+
  facet_wrap(~gene, ncol = 10)+
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("")+
  theme_bw()+
  ggtitle("Female Ovaries: Upregulated genes in NMX")+
  theme(legend.position = "none", text =element_text(size = 10))
p1


ggsave("imgs/FemaleOVT_Upreg.png", width = 9, height = 2.25)


x %>% 
  filter(gene %in% down_genes) -> xex
colnames(xex)

xex2 <- xex %>% pivot_longer(cols = 2:24, names_to = "ids")

p <- xex2 %>% full_join(a_data)

p$Group <- ifelse(p$Group == "Orange", "Control", "NMX")
p$group <- factor(p$Group, levels = c("Control", "NMX"))
p <- p %>% na.omit(.)

source('functions/geom_boxjitter.R')
library(viridis)

p2 <- ggplot(p, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA, alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE, position = position_dodge(0.85)) +
  scale_color_manual(values = c( "darkorange", "blue"))+
  scale_fill_manual(values = c( "darkorange", "blue"))+
  facet_wrap(~gene, ncol = 10)+
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("")+
  theme_bw()+
  ggtitle("Female Ovaries: Downregulated genes in NMX")+
  theme(legend.position = "none", text =element_text(size = 10))
p2


ggsave("imgs/FemaleVOVT_Downreg.png", width = 9, height = 2.25)

# boxplots for arc glu genes  genes 

ex <- readRDS("males/clean_data/limma_vdl_ARC.RDS")
head(ex)

x <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('gene') 


# sample information 
id <- read_csv("raw_data/traitData.csv")
colnames(id)
idx <- id %>% select(1:4,8:11)

#coldata arc
a_data <- idx
a_data$SampleName <- paste("A", a_data$SampleName, sep = "_")
a_data$ids <- a_data$SampleName
a_datax <- a_data %>% column_to_rownames(., var = "SampleName")
a_datax$bi_sex <- ifelse(a_datax$Sex == "M", 0, 1)
a_datax$bi_group <- ifelse(a_datax$Group == "Blue", 0, 1)


a_datax$Group



top_genes <- c("Bdnf", "Grin2b", "Grm1", "Shank1", "Gabra1", "Gabra2", "Gabra4")
down_genes <- c("Ghrh", "Lepr", "Gal", "Igfbp6", "Ttr", "Timp2", "Pomc", "Igfbp3", "Npy", "Thra","Thrb", "Th")

x %>% 
  filter(gene %in% down_genes) -> xex
colnames(xex)

xex2 <- xex %>% pivot_longer(cols = 2:14, names_to = "ids")

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
  ggtitle("Male ARC: Upregulated genes in NMX")+
  theme(legend.position = "none", text =element_text(size = 10))
p1


ggsave("imgs/MaleARC_Upreg.png", width = 10, height = 2.25)


p
ggplot(p, aes(P84,value, color = Group))+
  geom_point( size = 1) +
  stat_smooth(method = "lm", se = F)+
  facet_wrap(~gene, ncol = 6) +
  scale_y_continuous(expand = c(0, 1))+
  scale_color_manual(values = c( "darkorange", "blue"))+
  ylab("Normalized Expression") +
  xlab("P84 Body Weight (g)")+
  theme_minimal()+
  ggtitle("Male ARC")+
  theme( text =element_text(size = 15))


ggplot(p, aes(P1,value, color = Group))+
  geom_point( size = 1) +
  stat_smooth(method = "lm", se = F)+
  facet_wrap(~gene, ncol = 6) +
  scale_y_continuous(expand = c(0, 1))+
  scale_color_manual(values = c( "darkorange", "blue"))+
  ylab("Normalized Expression") +
  xlab("P1 Body Weight (g)")+
  theme_minimal()+
  ggtitle("Male ARC")+
  theme( text =element_text(size = 15))

p$gene

#Gal
gal <- p %>% filter(gene == "Gal")
g1 <- lm(value~group+P1, data = gal)
summary(g1)

g84 <- lm(value~group+P84, data = gal)
summary(g84)

#ghrh
gh <- p %>% filter(gene == "Ghrh")
gh1 <- lm(value~group+P1, data = gh)
summary(gh1) ##**

gh84 <- lm(value~group+P84, data = gh)

summary(gh84) ##**
cor.test(gh$P1,gh$value)
cor.test(gh$P84,gh$value)


#Igfbp3

ig3 <- p %>% filter(gene == "Igfbp3")
ig1 <- lm(value~group+P1, data = ig3)
summary(ig1)

ig84 <- lm(value~group+P84, data = ig3)
summary(ig84)  #0.0758

cor.test(ig3$P84,ig3$value) # 0.0597

ig84_inx <- lm(value~group*P84, data = ig3)
summary(ig84_inx)  #0.0758

#Igfbp6
ig6 <- p %>% filter(gene == "Igfbp6")
ig1x <- lm(value~group+P1, data = ig6)
summary(ig1x)

ig84x <- lm(value~groupP*84, data = ig6)
summary(ig84x)  #0.0740

cor.test(ig6$P1,ig6$value)

#Lepr
l <- p %>% filter(gene == "Lepr")
l1 <- lm(value~group+P1, data = l)
summary(ig1)

l84 <- lm(value~group+P84, data = l)
summary(l84) 


#Npy

n <- p %>% filter(gene == "Npy")
n1 <- lm(value~group+P1, data = n)
summary(n1)

n84 <- lm(value~group+P84, data = n)
summary(n84) 

n84_inx <- lm(value~group*P84, data = n)
summary(n84_inx) 

#Pomc

pp <- p %>% filter(gene == "Pomc")
p1 <- lm(value~group+P1, data = pp)
summary(p1)

p84 <- lm(value~group+P84, data = pp)
summary(p84)  #0.0758

#Th
th<- p %>% filter(gene == "Th")
th1 <- lm(value~group+P1, data = th)
summary(th1)

th84 <- lm(value~group+P84, data = th)
summary(th84)  #0.098


#Thra
#Thrb
#Timp2
#Ttr


tt <- p %>% filter(gene == "Ttr")
tt1 <- lm(value~group+P1, data = tt)
summary(tt1)

tt84 <- lm(value~group+P84, data = tt)
summary(tt84)  #0.0996

cor.test(tt$P84,tt$value)



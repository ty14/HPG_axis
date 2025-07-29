library(tidyverse)


df <- read_csv("males/clean_data/Top_GOTerms_males.csv")
head(df)


a.df <- df %>% filter(tissue == "ARC")
colnames(a.df)
a.df$eFDR <- a.df$qvalue

m_a <- ggplot(a.df, aes(x = GeneRatio, y = Description, size = -log10(eFDR))) +
  geom_point(aes(fill = direction), shape = 21,
    color = "black",
    stroke = 0.8, show.legend = T) +
  labs(title = "A. Male ARC", x = "Gene Ratio", y = "") +
  theme_bw()  +
  labs(color = "Regulation") +
  scale_fill_viridis_d(name = "Regulation", option = "E", direction = -1)+
  scale_size_continuous(limits = c(0,25), 
                        breaks = c(5, 7, 10))+
  scale_x_continuous(
    limits = c(0.01,0.11),
    breaks = c(0.02,0.04,0.06,0.08,0.10))+
  theme(
    # axis.line = element_line(colour = 'black', size = 1),
    axis.ticks = element_line(colour = "black", size = 1),
    text = element_text(size=20, colour="black"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  ) 

ggsave("imgs/maleARC_go_bubble_sex.png", m_a, width = 8, height = 4, dpi = 600)

v.df <- df %>% filter(tissue == "VMN")
colnames(v.df)
v.df$eFDR <- v.df$qvalue

m_v <- ggplot(v.df, aes(x = GeneRatio, y = Description, size = -log10(eFDR))) +
  geom_point(aes(fill = direction), shape = 21,
             color = "black",
             stroke = 0.8, show.legend = T)  +
  labs(title = "C. Male VMN", x = "Gene Ratio", y = "") +
  theme_bw()  +
  labs(color = "Regulation") +
  scale_fill_viridis_d(name = "Regulation", option = "E", direction = -1)+
  scale_size_continuous(limits = c(0,25), 
                        breaks = c(5, 7, 10))+
  scale_x_continuous(
    limits = c(0.01,0.11),
    breaks = c(0.02,0.04,0.06,0.08,0.10))+
  theme(
    # axis.line = element_line(colour = 'black', size = 1),
    axis.ticks = element_line(colour = "black", size = 1),
    text = element_text(size=20, colour="black"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  ) 

m_v
ggsave("imgs/maleVMN_go_bubble_use.png", m_v, width = 8, height = 4, dpi = 600)


p.df <- df %>% filter(tissue == "PIT")
colnames(p.df)
p.df$eFDR <- p.df$qvalue

m_p <- ggplot(p.df, aes(x = GeneRatio, y = Description, size = -log10(eFDR))) +
  geom_point(aes(fill = direction), shape = 21,
             color = "black",
             stroke = 0.8, show.legend = T) +
  labs(title = "E. Male PIT", x = "Gene Ratio", y = "") +
  theme_bw()  +
  labs(color = "Regulation") +
  scale_fill_viridis_d(name = "Regulation", option = "E", direction = 1)+
  scale_size_continuous(limits = c(0,15), 
                        breaks = c(5, 7, 10))+
  scale_x_continuous(
    limits = c(0.01,0.11),
    breaks = c(0.02,0.04,0.06,0.08,0.10))+
  theme(
    # axis.line = element_line(colour = 'black', size = 1),
    axis.ticks = element_line(colour = "black", size = 1),
    text = element_text(size=20, colour="black"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  ) 


m_p
ggsave("imgs/malePIT_go_bubble_use.png", m_p, width = 8, height = 4, dpi = 600)


### female hypo and pit 

df <- read_csv("females/clean_data/Top_GOTerms_females.csv")
head(df)


a.df <- df %>% filter(tissue == "ARC")
colnames(a.df)
a.df$eFDR <- a.df$qvalue

fe.a <- ggplot(a.df, aes(x = GeneRatio, y = Description, size = -log10(eFDR))) +
  geom_point(aes(fill = direction), shape = 21,
             color = "black",
             stroke = 0.8, show.legend = T) +
  labs(title = "B. Female ARC", x = "Gene Ratio", y = "") +
  theme_bw()  +
  labs(color = "Regulation") +
  scale_fill_viridis_d(name = "Regulation", option = "E", direction = -1)+
  scale_size_continuous(limits = c(0,25), 
                        breaks = c(5, 7, 10))+
  scale_x_continuous(
    limits = c(0.01,0.11),
    breaks = c(0.02,0.04,0.06,0.08,0.10)) +
  theme(
    # axis.line = element_line(colour = 'black', size = 1),
    axis.ticks = element_line(colour = "black", size = 1),
    text = element_text(size=20, colour="black"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  ) 


ggsave("imgs/femaleARC_go_bubble_use.png", fe.a, width =7.35, height = 4, dpi = 600)


v.df <- df %>% filter(tissue == "VMN")
colnames(v.df)
v.df$eFDR <- v.df$qvalue

fe.v <- ggplot(v.df, aes(x = GeneRatio, y = Description, size = -log10(eFDR))) +
  geom_point(aes(fill = direction), shape = 21,
             color = "black",
             stroke = 0.8, show.legend = T) +
  labs(title = "D. Female VMN", x = "Gene Ratio", y = "") +
  theme_bw()  +
  labs(color = "Regulation") +
  scale_fill_viridis_d(name = "Regulation", option = "E", direction = -1)+
  scale_size_continuous(limits = c(0,25), 
                        breaks = c(5, 7, 10))+
  scale_x_continuous(
    limits = c(0.01,0.11),
    breaks = c(0.02,0.04,0.06,0.08,0.1))+
  theme(
    
    # axis.line = element_line(colour = 'black', size = 1),
    axis.ticks = element_line(colour = "black", size = 1),
    text = element_text(size=20, colour="black"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  ) 


ggsave("imgs/femaleVMN_go_bubble_use.png", fe.v, width = 7.2, height = 4, dpi = 600)


p.df <- df %>% filter(tissue == "PIT")
colnames(p.df)
p.df$eFDR <- p.df$qvalue

fe.p <- ggplot(p.df, aes(x = GeneRatio, y = Description, size = -log10(eFDR))) +
  geom_point(aes(fill = direction), shape = 21,
             color = "black",
             stroke = 0.8, show.legend = T) +
  labs(title = "F. Female PIT", x = "Gene Ratio", y = "") +
  theme_bw()  +
  labs(color = "Regulation") +
  scale_fill_viridis_d(name = "Regulation", option = "E", direction = -1)+
  scale_size_continuous(limits = c(0,15), 
                        breaks = c(5, 7, 10))+
  scale_x_continuous(
    limits = c(0.01,0.11),
    breaks = c(0.02,0.04,0.06,0.08,0.10))+
  theme(
    # axis.line = element_line(colour = 'black', size = 1),
    axis.ticks = element_line(colour = "black", size = 1),
    text = element_text(size=20, colour="black"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  )  



ggsave("imgs/femalePIT_go_bubble_use.png", fe.p, width = 7.75, height = 4, dpi = 600)




x <- gridExtra::grid.arrange(m_a, m_v, m_p, fe.a, fe.v,fe.p, ncol =3)

ggsave("imgs/go_bubbleplots.png", x, width = 20, height = 10, dpi = 600)






# ##################### GONADS 

df <- read_csv("males/clean_data/Top_GOTerms_males.csv")
head(df)


t.df <- df %>% filter(tissue == "TES")
colnames(p.df)
t.df$eFDR <- t.df$qvalue

tp <- ggplot(t.df, aes(x = GeneRatio, y = Description, size = -log10(eFDR))) +
  geom_point(aes(fill = direction), shape = 21,
             color = "black",
             stroke = 0.8, show.legend = T) +
  labs(title = "A. Testes", x = "Gene Ratio", y = "") +
  theme_bw()  +
  labs(color = "Regulation") +
  scale_fill_viridis_d(name = "Regulation", option = "E", direction = -1)+
  scale_size_continuous(limits = c(0,15), 
                        breaks = c(5, 7, 10))+
  scale_x_continuous(
    limits = c(0.01,0.11),
    breaks = c(0.02,0.04,0.06,0.08,0.10))+
  theme(
    # axis.line = element_line(colour = 'black', size = 1),
    axis.ticks = element_line(colour = "black", size = 1),
    text = element_text(size=20, colour="black"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  )  

ggsave("imgs/testes_go_bubble_use.png", tp, width = 6.4, height =3.5 , dpi = 600)

df <- read_csv("females/clean_data/Top_GOTerms_females.csv")
head(df)


t.df <- df %>% filter(tissue == "OVT")
colnames(p.df)
t.df$eFDR <- t.df$qvalue

ovt <- ggplot(t.df, aes(x = GeneRatio, y = Description, size = -log10(eFDR))) +
  geom_point(aes(fill = direction), shape = 21,
             color = "black",
             stroke = 0.8, show.legend = T) +
  labs(title = "B. Ovaries", x = "Gene Ratio", y = "") +
  theme_bw()  +
  labs(color = "Regulation") +
  scale_fill_viridis_d(name = "Regulation", option = "E", direction = -1)+
  scale_size_continuous(limits = c(0,15), 
                        breaks = c(5, 7, 10))+
  scale_x_continuous(
    limits = c(0.01,0.11),
    breaks = c(0.02,0.04,0.06,0.08,0.10))+
  theme(
    # axis.line = element_line(colour = 'black', size = 1),
    axis.ticks = element_line(colour = "black", size = 1),
    text = element_text(size=20, colour="black"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  )  

ggsave("imgs/ovaries_go_bubble.png", ovt, width = 7, height = 3.5, dpi = 600)


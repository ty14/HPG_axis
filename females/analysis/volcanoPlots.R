#female VPlots
library(EnhancedVolcano)
library(tidyverse)

##ARC
fe_a  <- readRDS("males/clean_data/limma_eFDR_ARC.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val) %>% distinct(.) 

head(fe_a)
dc <- fe_a %>% mutate(contrast = "NMX vs Control") %>% mutate(log10 = -log10(P.Value))

dc$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
dc$diffexpressed[dc$logFC > 0.2 & dc$P.Value < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dc$diffexpressed[dc$logFC < -0.2 & dc$P.Value < 0.05] <- "DOWN"


dcx <- dc %>% filter(.,logFC >= 1.1)%>% filter(P.Value < 0.05)%>% filter(!grepl("LOC", genes))
dcxx <- dc %>% filter(.,logFC <= -1.44) %>% filter(P.Value < 0.05) %>% filter(!grepl("LOC", genes))
dc$log10 <- ifelse(dc$log10 == Inf, 4,dc$log10)
# dc$log10 <- ifelse(dc$log10 == Inf, 4,dc$loggenes10)
# dcx <- dc %>% filter(.,abs(logFC) >= 1.5) %>% filter(P.Value < 0.05) 

vp_csub <- ggplot(data = dc, 
                  aes(x = logFC, 
                      y = log10, 
                      colour=diffexpressed)) +
  geom_point(alpha=0.25, size=3.5) +
  scale_color_manual(values=c("orange", "grey","blue"))+
  xlim(c(-2.5, 2.5)) +
  ylim(0,4)+
  # geom_vline(xintercept=c(-0.2,0.2),lty=4,col="black",lwd=0.8) +
  # geom_vline(xintercept=c(-.75,.75),lty=4,col="black",lwd=0.8) +
  # geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = genes), color = "black",  hjust = -0.15,vjust =0.9, size = 5)+
  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = genes), color = "black", hjust = 1, vjust = -1, size = 5)+
  labs(x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR)))+
  theme_bw() +
  annotate(geom="text", x=2, y=.5, label=paste0("Up ", "\n", "in NMX"),
           color="black", size = 8)+
  annotate(geom="text", x=-2, y=.5, label=paste0("Down ", "\n", "in NMX"),
           color="black",, size = 8)+
  scale_x_continuous(limits = c(-3,3),breaks = c(-2.5,-1.5,-0.5,0,0.5,1.5,2.5),
                     labels = c("-2.5", "-1.5", "-0.5", "0", "0.5", "1.5", "2.5"))+
  # ggtitle("ARC")+
  theme(axis.text.x = element_text(vjust = 1,size = 25),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 25),
        axis.text = element_text(color="#3C3C3C",size = 25),
        axis.title = element_text(size = 25),   
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 30),
        # text = element_text(size = 25)
        
  )

vp_csub

ggsave("imgs/male_arc_vplot.png", vp_csub, height = 6, width = 6.5, dpi = 300)


###VMN
fe_v  <- readRDS("males/clean_data/limma_eFDR_VMN.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val) %>% distinct(.) 

head(fe_v)
dc <- fe_v %>% mutate(contrast = "NMX vs Control") %>% mutate(log10 = -log10(P.Value))

dc$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
dc$diffexpressed[dc$logFC > 0.2 & dc$P.Value < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dc$diffexpressed[dc$logFC < -0.2 & dc$P.Value < 0.05] <- "DOWN"


dcx <- dc %>% filter(.,logFC >= 1.5)%>% filter(P.Value < 0.05)
dcxx <- dc %>% filter(.,logFC <= -1.1) %>% filter(P.Value < 0.05)
dc$log10 <- ifelse(dc$log10 == Inf, 4,dc$log10)
# dcx <- dc %>% filter(.,abs(logFC) >= 1.5) %>% filter(P.Value < 0.05) 

vp_csub <- ggplot(data = dc, 
                  aes(x = logFC, 
                      y = log10, 
                      colour=diffexpressed)) +
  geom_point(alpha=0.25, size=3.5) +
  scale_color_manual(values=c("orange", "grey","blue"))+
  xlim(c(-2.5, 2.5)) +
  ylim(0,4)+
  # geom_vline(xintercept=c(-0.2,0.2),lty=4,col="black",lwd=0.8) +
  # geom_vline(xintercept=c(-.75,.75),lty=4,col="black",lwd=0.8) +
  # geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = genes), color = "black",  hjust = -0.5,vjust =-0.8, size = 5)+
  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = genes), color = "black", hjust = 0.8, vjust =0.8, size = 5)+
  labs(x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR)))+
  theme_bw() +
  annotate(geom="text", x=2, y=.5, label=paste0("Up ", "\n", "in NMX"),
           color="black", size = 8)+
  annotate(geom="text", x=-2, y=.5, label=paste0("Down ", "\n", "in NMX"),
           color="black",, size = 8)+
  scale_x_continuous(limits = c(-3,3),breaks = c(-2.5,-1.5,-0.5,0,0.5,1.5,2.5),
                     labels = c("-2.5", "-1.5", "-0.5", "0", "0.5", "1.5", "2.5"))+
  # ggtitle("VMN")+
  theme(axis.text.x = element_text(vjust = 1,size = 25),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 25),
        axis.text = element_text(color="#3C3C3C",size = 25),
        axis.title = element_text(size = 25),   
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 30),
        # text = element_text(size = 25)
        
  )

vp_csub

ggsave("imgs/male_vmn_vplot.png", vp_csub, height = 6, width = 6.5, dpi = 300)

###VMN
fe_v  <- readRDS("males/clean_data/limma_eFDR_PIT.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val) %>% distinct(.) 

head(fe_v)
dc <- fe_v %>% mutate(contrast = "NMX vs Control") %>% mutate(log10 = -log10(P.Value))

dc$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
dc$diffexpressed[dc$logFC > 0.2 & dc$P.Value < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dc$diffexpressed[dc$logFC < -0.2 & dc$P.Value < 0.05] <- "DOWN"


dcx <- dc %>% filter(.,logFC >= 0.6)%>% filter(P.Value < 0.05)
dcxx <- dc %>% filter(.,logFC <= -0.65) %>% filter(P.Value < 0.05)
dc$log10 <- ifelse(dc$log10 == Inf, 4,dc$log10)
# dcx <- dc %>% filter(.,abs(logFC) >= 1.5) %>% filter(P.Value < 0.05) 

vp_csub <- ggplot(data = dc, 
                  aes(x = logFC, 
                      y = log10, 
                      colour=diffexpressed)) +
  geom_point(alpha=0.25, size=3.5) +
  scale_color_manual(values=c("orange", "grey","blue"))+
  xlim(c(-2.5, 2.5)) +
  ylim(0,4)+
  # geom_vline(xintercept=c(-0.2,0.2),lty=4,col="black",lwd=0.8) +
  # geom_vline(xintercept=c(-.75,.75),lty=4,col="black",lwd=0.8) +
  # geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = genes), color = "black",  hjust = -0.75,vjust =-0.8, size = 5)+
  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = genes), color = "black", hjust = 0.8, vjust =0.8, size = 5)+
  labs(x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR)))+
  theme_bw() +
  annotate(geom="text", x=2, y=.5, label=paste0("Up ", "\n", "in NMX"),
           color="black", size = 8)+
  annotate(geom="text", x=-2, y=.5, label=paste0("Down ", "\n", "in NMX"),
           color="black",, size = 8)+
  scale_x_continuous(limits = c(-3,3),breaks = c(-2.5,-1.5,-0.5,0,0.5,1.5,2.5),
                     labels = c("-2.5", "-1.5", "-0.5", "0", "0.5", "1.5", "2.5"))+
  # ggtitle("PIT")+
  theme(axis.text.x = element_text(vjust = 1,size = 25),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 25),
        axis.text = element_text(color="#3C3C3C",size = 25),
        axis.title = element_text(size = 25),   
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 30),
        # text = element_text(size = 25)
        
  )

vp_csub

ggsave("imgs/male_pit_vplot.png", vp_csub, height = 6, width = 6.5, dpi = 300)


###VMN
fe_v  <- readRDS("males/clean_data/limma_eFDR_Testis.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val) %>% distinct(.) 

head(fe_v)
dc <- fe_v %>% mutate(contrast = "NMX vs Control") %>% mutate(log10 = -log10(P.Value))

dc$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
dc$diffexpressed[dc$logFC > 0.2 & dc$P.Value < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dc$diffexpressed[dc$logFC < -0.2 & dc$P.Value < 0.05] <- "DOWN"


dcx <- dc %>% filter(.,logFC >= 0.4)%>% filter(P.Value < 0.05) %>% filter(!grepl("LOC", genes))
dcxx <- dc %>% filter(.,logFC <= -0.68) %>% filter(P.Value < 0.05)
dc$log10 <- ifelse(dc$log10 == Inf, 4,dc$log10)
# dcx <- dc %>% filter(.,abs(logFC) >= 1.5) %>% filter(P.Value < 0.05) 

vp_csub <- ggplot(data = dc, 
                  aes(x = logFC, 
                      y = log10, 
                      colour=diffexpressed)) +
  geom_point(alpha=0.25, size=3.5) +
  scale_color_manual(values=c("orange", "grey","blue"))+
  xlim(c(-2.5, 2.5)) +
  ylim(0,4)+
  # geom_vline(xintercept=c(-0.2,0.2),lty=4,col="black",lwd=0.8) +
  # geom_vline(xintercept=c(-.75,.75),lty=4,col="black",lwd=0.8) +
  # geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = genes), color = "black",  hjust = -0.75,vjust =-0.8, size = 5)+
  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = genes), color = "black", hjust = 0.8, vjust =0.8, size = 5)+
  labs(x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR)))+
  theme_bw() +
  annotate(geom="text", x=2, y=.5, label=paste0("Up ", "\n", "in NMX"),
           color="black", size = 8)+
  annotate(geom="text", x=-2, y=.5, label=paste0("Down ", "\n", "in NMX"),
           color="black",, size = 8)+
  scale_x_continuous(limits = c(-3,3),breaks = c(-2.5,-1.5,-0.5,0,0.5,1.5,2.5),
                     labels = c("-2.5", "-1.5", "-0.5", "0", "0.5", "1.5", "2.5"))+
  # ggtitle("Testes")+
  theme(axis.text.x = element_text(vjust = 1,size = 25),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 25),
        axis.text = element_text(color="#3C3C3C",size = 25),
        axis.title = element_text(size = 25),   
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 25),
        # text = element_text(size = 25)
        
  )

vp_csub

ggsave("imgs/male_tes_vplot_endo.png", vp_csub, height = 6, width = 6.5, dpi = 300)

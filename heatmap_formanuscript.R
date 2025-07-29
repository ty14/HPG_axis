library(tidyverse)

###get DEGs
### get all the genes for each tissue type:
#females 
fe_vmn <- readRDS("females/clean_data/limma_eFDR_VMN.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val) %>% distinct(.)

fe_arc <- readRDS("females/clean_data/limma_eFDR_ARC.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val) %>% distinct(.)

fe_pit <- readRDS("females/clean_data/limma_eFDR_PIT.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val) %>% distinct(.)

fe_ovt <- readRDS("females/clean_data/limma_eFDR_OVT.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val) %>% distinct(.)


m_vmn <- readRDS("males/clean_data/limma_eFDR_VMN.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val) %>% distinct(.)

m_arc <- readRDS("males/clean_data/limma_eFDR_ARC.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val) %>% distinct(.)

m_pit <- readRDS("males/clean_data/limma_eFDR_PIT.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val) %>% distinct(.)

m_tes <- readRDS("males/clean_data/limma_eFDR_Testis.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val) %>% distinct(.)


#endocrine genes 
#Esr1, Ar, Gal, Esr2, Cyp19a1, Egfr, Igf1,   Kiss1, Pdyn, Tac3, Lhb, Fshb, Nr5a2, Hmgcr, Igf1r

#synaptic plasticity 
#Bdfn, Grm1, Grm4, Gabra1, Gabra2, Gabra4, Synpo1, Oxt, Shank1, Grip1

#metabolic genes 
#Ghrh, Lepr, Igfbp6, Igfbp3, Ttr, Cck, Pmch, Pomc,Igf1r, Timp1, Npy, Agrp, Hcrt1, Hcrt2,  Tshm Tshr, Thra, Thrb, Trh,  

#Cilia genes 
# Celsr3, Cfap65, Dnah11, Dnah12, Dnai7, Erich3, Dnah10, Cfap161, Spef1, Ccdc39, Ccdc40

#immune genes 
#Il6r, Camk2d, Sox2, RT1-Da, RT1-Db1, Ccl5,Cc17, RT1-DMa, Rt1-DMb, Mpeg, Ptgds Tnf, Ifna, Nfkb, Il6, Il1, Fam107a, Wnt4, Wnt 6, Il10, Nfkbia, Ptgds

#epigenetics/RNA processing 
# Kmt5a , Dnmt3a, Dnmt1,Hdac1, Hdac10, Hdac7, Trim65,Trim45, Trim28

#Glucocorticoids
# Fkbp5, Nr3c1, Hsp90AA1, Hsp70, P23, Ptges2, G6pc, Foxo1, Dusp1, Crh, Avp, Hsb11b1, Hsb11b2, Nts



#first estrogen genes

es <- c("Esr1", "Ar", "Gal", "Esr2", "Cyp17a1", "Cyp19a1", "Egfr", "Igf1","Igf1r", "Kiss1", "Pdyn", "Tac3", "Lhb", "Fshb", "Nr5a2", "Hmgcr")

fe.arc.es <- fe_arc %>% filter(genes %in% es) %>% select(genes, "F-ARC" = logFC) %>% unique(.) %>% arrange(-`F-ARC`)
m.arc.es <- m_arc %>% filter(genes %in% es)%>% select(genes, "M-ARC" = logFC) %>% unique(.) %>% arrange(-`M-ARC`)

fe.vmn.es <- fe_vmn %>% filter(genes %in% es)%>% select(genes, "F-VMN" = logFC) %>% unique(.) %>% arrange(-`F-VMN`)
m.vmn.es <- m_vmn %>% filter(genes %in% es)%>% select(genes, "M-VMN" = logFC) %>% unique(.) %>% arrange(-`M-VMN`)

fe.pit.es <- fe_pit %>% filter(genes %in% es)%>% select(genes, "F-PIT" = logFC) %>% unique(.) %>% arrange(-`F-PIT`)
m.pit.es <- m_pit %>% filter(genes %in% es)%>% select(genes, "M-PIT" = logFC) %>% unique(.) %>% arrange(-`M-PIT`)

ovt.es <- fe_ovt %>% filter(genes %in% es)%>% select(genes, "Ovaries" = logFC) %>% unique(.) %>% arrange(-Ovaries)
tes.es <- m_tes %>% filter(genes %in% es)%>% select(genes, "Testes" = logFC) %>% unique(.) %>% arrange(-`Testes`)

hm_esr <-m.arc.es %>% full_join(fe.arc.es)  %>% full_join(m.vmn.es)%>% full_join(fe.vmn.es) %>% full_join(m.pit.es) %>% full_join(fe.pit.es)  %>% full_join(tes.es)%>% full_join(ovt.es) 
colnames(hm_esr)
hm_esrx <- hm_esr %>% pivot_longer(cols = 2:9, names_to = "tissue", values_to = "NXM vs. Control")
hm_esrx[is.na(hm_esrx)] <- 0
hm_esrx


library(ComplexHeatmap)
hm_esrxx <- hm_esr%>% arrange(genes) %>% column_to_rownames(., var= "genes") %>% as.matrix()
hm_esrxx[is.na(hm_esrxx)] <- 0

Heatmap(hm_esrxx, name = "LogFC",cluster_columns = FALSE,cluster_rows = FALSE,
        show_heatmap_legend = T, column_names_rot = 60, rect_gp = gpar(col = "white", lwd = 2))

#metabolic genes 

m <- c("Ghrh", "Lepr", "Igfbp6", "Igfbp3", "Ttr", "Cck", "Pmch", "Pomc", "Timp1", "Npy", "Thra", "Thrb") 

fe.arc.m <- fe_arc %>% filter(genes %in% m) %>% select(genes, "F-ARC" = logFC) %>% unique(.) %>% arrange(-`F-ARC`)
m.arc.m <- m_arc %>% filter(genes %in% m)%>% select(genes, "M-ARC" = logFC) %>% unique(.) %>% arrange(-`M-ARC`)

fe.vmn.m <- fe_vmn %>% filter(genes %in% m)%>% select(genes, "F-VMN" = logFC) %>% unique(.) %>% arrange(-`F-VMN`)
m.vmn.m <- m_vmn %>% filter(genes %in% m)%>% select(genes, "M-VMN" = logFC) %>% unique(.) %>% arrange(-`M-VMN`)

fe.pit.m <- fe_pit %>% filter(genes %in% m)%>% select(genes, "F-PIT" = logFC) %>% unique(.) %>% arrange(-`F-PIT`)
m.pit.m <- m_pit %>% filter(genes %in% m)%>% select(genes, "M-PIT" = logFC) %>% unique(.) %>% arrange(-`M-PIT`)

ovt.m <- fe_ovt %>% filter(genes %in% m)%>% select(genes, "Ovaries" = logFC) %>% unique(.) %>% arrange(-Ovaries)
tes.m <- m_tes %>% filter(genes %in% m)%>% select(genes, "Testes" = logFC) %>% unique(.) %>% arrange(-`Testes`)

hm_m <-m.arc.m %>% full_join(fe.arc.m)  %>% full_join(m.vmn.m)%>% full_join(fe.vmn.m) %>% full_join(m.pit.m) %>% full_join(fe.pit.m)  %>% full_join(tes.m)%>% full_join(ovt.m) 
colnames(hm_m)
hm_mx <- hm_m %>% pivot_longer(cols = 2:9, names_to = "tissue", values_to = "NXM vs. Control")
hm_mx[is.na(hm_mx)] <- 0
hm_mx


library(ComplexHeatmap)
hm_mxx <- hm_m%>% arrange(genes) %>% column_to_rownames(., var= "genes") %>% as.matrix()
hm_mxx[is.na(hm_mxx)] <- 0

Heatmap(hm_mxx, name = "LogFC",cluster_columns = FALSE,cluster_rows = FALSE,
        show_heatmap_legend = T, column_names_rot = 60, rect_gp = gpar(col = "white", lwd = 2))


#Synpatic genes  

s <- c("Bdnf", "Grm1", "Grm4", "Gabra1", "Gabra2", "Gabra4", "Synpo1", "Oxt", "Shank1", "Grip1") 

fe.arc.s <- fe_arc %>% filter(genes %in% s) %>% select(genes, "F-ARC" = logFC) %>% unique(.) %>% arrange(-`F-ARC`)
m.arc.s <- m_arc %>% filter(genes %in% s)%>% select(genes, "M-ARC" = logFC) %>% unique(.) %>% arrange(-`M-ARC`)

fe.vmn.s <- fe_vmn %>% filter(genes %in% s)%>% select(genes, "F-VMN" = logFC) %>% unique(.) %>% arrange(-`F-VMN`)
m.vmn.s <- m_vmn %>% filter(genes %in% s)%>% select(genes, "M-VMN" = logFC) %>% unique(.) %>% arrange(-`M-VMN`)

fe.pit.s <- fe_pit %>% filter(genes %in% s)%>% select(genes, "F-PIT" = logFC) %>% unique(.) %>% arrange(-`F-PIT`)
m.pit.s <- m_pit %>% filter(genes %in% s)%>% select(genes, "M-PIT" = logFC) %>% unique(.) %>% arrange(-`M-PIT`)

ovt.s <- fe_ovt %>% filter(genes %in% s)%>% select(genes, "Ovaries" = logFC) %>% unique(.) %>% arrange(-Ovaries)
tes.s <- m_tes %>% filter(genes %in% s)%>% select(genes, "Testes" = logFC) %>% unique(.) %>% arrange(-`Testes`)

hm_s <-m.arc.s %>% full_join(fe.arc.s)  %>% full_join(m.vmn.s)%>% full_join(fe.vmn.s) %>% full_join(m.pit.s) %>% full_join(fe.pit.s)  %>% full_join(tes.s)%>% full_join(ovt.s) 
colnames(hm_s)
hm_sx <- hm_s %>% pivot_longer(cols = 2:9, names_to = "tissue", values_to = "NXM vs. Control")
hm_sx[is.na(hm_sx)] <- 0
hm_sx


library(ComplexHeatmap)
hm_sxx <- hm_s%>% arrange(genes) %>% column_to_rownames(., var= "genes") %>% as.matrix()
hm_sxx[is.na(hm_sxx)] <- 0

Heatmap(hm_sxx, name = "LogFC",cluster_columns = FALSE,cluster_rows = FALSE,
        show_heatmap_legend = T, column_names_rot = 60, rect_gp = gpar(col = "white", lwd = 2))


#Cilia genes 
c <- c("Celsr3", "Cfap65", "Dnah11", "Dnah12", "Dnai7", "Erich3", "Dnah10", "Cfap161", "Spef1", "Ccdc39", "Ccdc40")

fe.arc.c <- fe_arc %>% filter(genes %in% c) %>% select(genes, "F-ARC" = logFC) %>% unique(.) %>% arrange(-`F-ARC`)
m.arc.c<- m_arc %>% filter(genes %in% c)%>% select(genes, "M-ARC" = logFC) %>% unique(.) %>% arrange(-`M-ARC`)

fe.vmn.c <- fe_vmn %>% filter(genes %in% c)%>% select(genes, "F-VMN" = logFC) %>% unique(.) %>% arrange(-`F-VMN`)
m.vmn.c <- m_vmn %>% filter(genes %in% c)%>% select(genes, "M-VMN" = logFC) %>% unique(.) %>% arrange(-`M-VMN`)

fe.pit.c <- fe_pit %>% filter(genes %in% c)%>% select(genes, "F-PIT" = logFC) %>% unique(.) %>% arrange(-`F-PIT`)
m.pit.c <- m_pit %>% filter(genes %in% c)%>% select(genes, "M-PIT" = logFC) %>% unique(.) %>% arrange(-`M-PIT`)

ovt.c <- fe_ovt %>% filter(genes %in% c)%>% select(genes, "Ovaries" = logFC) %>% unique(.) %>% arrange(-Ovaries)
tes.c <- m_tes %>% filter(genes %in% c)%>% select(genes, "Testes" = logFC) %>% unique(.) %>% arrange(-`Testes`)

hm_c <-m.arc.c %>% full_join(fe.arc.c)  %>% full_join(m.vmn.c)%>% full_join(fe.vmn.c) %>% full_join(m.pit.c) %>% full_join(fe.pit.c)  %>% full_join(tes.c)%>% full_join(ovt.c) 
colnames(hm_c)
hm_cx <- hm_c %>% pivot_longer(cols = 2:9, names_to = "tissue", values_to = "NXM vs. Control")
hm_cx[is.na(hm_cx)] <- 0
hm_cx


#immune genes 
im <- c("Il6r", "Camk2d", "Sox2", "RT1-Da", "RT1-Db1", "Ccl5","Cc17", "RT1-DMa", "Rt1-DMb", "Mpeg", "Tnfa", "Ifna", "Nfkb", "Il6", "Il1", "Fam107a", "Wnt4", "Wnt6", "Il10", "Nfkbia", "Ptgds", "Ptges2")

fe.arc.im <- fe_arc %>% filter(genes %in% im) %>% select(genes, "F-ARC" = logFC) %>% unique(.) %>% arrange(-`F-ARC`)
m.arc.im<- m_arc %>% filter(genes %in% im)%>% select(genes, "M-ARC" = logFC) %>% unique(.) %>% arrange(-`M-ARC`)

fe.vmn.im <- fe_vmn %>% filter(genes %in% im)%>% select(genes, "F-VMN" = logFC) %>% unique(.) %>% arrange(-`F-VMN`)
m.vmn.im <- m_vmn %>% filter(genes %in% im)%>% select(genes, "M-VMN" = logFC) %>% unique(.) %>% arrange(-`M-VMN`)

fe.pit.im <- fe_pit %>% filter(genes %in% im)%>% select(genes, "F-PIT" = logFC) %>% unique(.) %>% arrange(-`F-PIT`)
m.pit.im <- m_pit %>% filter(genes %in% im)%>% select(genes, "M-PIT" = logFC) %>% unique(.) %>% arrange(-`M-PIT`)

ovt.im <- fe_ovt %>% filter(genes %in% im)%>% select(genes, "Ovaries" = logFC) %>% unique(.) %>% arrange(-Ovaries)
tes.im <- m_tes %>% filter(genes %in% im)%>% select(genes, "Testes" = logFC) %>% unique(.) %>% arrange(-`Testes`)

hm_im <-m.arc.im %>% full_join(fe.arc.im)  %>% full_join(m.vmn.im)%>% full_join(fe.vmn.im) %>% full_join(m.pit.im) %>% full_join(fe.pit.im)  %>% full_join(tes.im)%>% full_join(ovt.im) 
colnames(hm_im)
hm_imx <- hm_im %>% pivot_longer(cols = 2:9, names_to = "tissue", values_to = "NXM vs. Control")
# hm_imx[is.na(
library(ComplexHeatmap)
hm_imxx <- hm_im%>% arrange(genes) %>% column_to_rownames(., var= "genes") %>% as.matrix()
hm_imxx[is.na(hm_imxx)] <- 0

Heatmap(hm_imxx, name = "LogFC",cluster_columns = FALSE,cluster_rows = FALSE,
        show_heatmap_legend = T, column_names_rot = 60, rect_gp = gpar(col = "white", lwd = 2))



#epigenetics/RNA processing 
ep <- c("Kmt5a", "Dnmt3a", "Dnmt1","Hdac1", "Hdac10", "Hdac7", "Trim65","Trim45", "Trim28", "Kmt1a", "Kmt2a", "Kmt3a", "Kmt2b")

fe.arc.ep <- fe_arc %>% filter(genes %in% ep) %>% select(genes, "F-ARC" = logFC) %>% unique(.) %>% arrange(-`F-ARC`)
m.arc.ep<- m_arc %>% filter(genes %in% ep)%>% select(genes, "M-ARC" = logFC) %>% unique(.) %>% arrange(-`M-ARC`)

fe.vmn.ep <- fe_vmn %>% filter(genes %in% ep)%>% select(genes, "F-VMN" = logFC) %>% unique(.) %>% arrange(-`F-VMN`)
m.vmn.ep <- m_vmn %>% filter(genes %in% ep)%>% select(genes, "M-VMN" = logFC) %>% unique(.) %>% arrange(-`M-VMN`)

fe.pit.ep <- fe_pit %>% filter(genes %in% ep)%>% select(genes, "F-PIT" = logFC) %>% unique(.) %>% arrange(-`F-PIT`)
m.pit.ep <- m_pit %>% filter(genes %in% ep)%>% select(genes, "M-PIT" = logFC) %>% unique(.) %>% arrange(-`M-PIT`)

ovt.ep <- fe_ovt %>% filter(genes %in% ep)%>% select(genes, "Ovaries" = logFC) %>% unique(.) %>% arrange(-Ovaries)
tes.ep <- m_tes %>% filter(genes %in% ep)%>% select(genes, "Testes" = logFC) %>% unique(.) %>% arrange(-`Testes`)

hm_ep <-m.arc.ep %>% full_join(fe.arc.ep)  %>% full_join(m.vmn.ep)%>% full_join(fe.vmn.ep) %>% full_join(m.pit.ep) %>% full_join(fe.pit.ep)  %>% full_join(tes.ep)%>% full_join(ovt.ep) 
colnames(hm_ep)
hm_epx <- hm_ep %>% pivot_longer(cols = 2:9, names_to = "tissue", values_to = "NXM vs. Control")
# hm_epx[is.na(
library(ComplexHeatmap)
hm_epxx <- hm_ep%>% arrange(genes) %>% column_to_rownames(., var= "genes") %>% as.matrix()
hm_epxx[is.na(hm_epxx)] <- 0

Heatmap(hm_epxx, name = "LogFC",cluster_columns = FALSE,cluster_rows = FALSE,
        show_heatmap_legend = T, column_names_rot = 60, rect_gp = gpar(col = "white", lwd = 2))




#gluc
g <- c("Fkbp5", "Nr3c1", "Hsp90AA1", "Hsp70","Ch", "G6pc", "Foxo1", "Dusp1", "Crh", "Avp","Crhr", "Hsb11b1", "Hsb11b2", "Nts")
fe.arc.g <- fe_arc %>% filter(genes %in% g) %>% select(genes, "F-ARC" = logFC) %>% unique(.) %>% arrange(-`F-ARC`)
m.arc.g<- m_arc %>% filter(genes %in% g)%>% select(genes, "M-ARC" = logFC) %>% unique(.) %>% arrange(-`M-ARC`)

fe.vmn.g <- fe_vmn %>% filter(genes %in% g)%>% select(genes, "F-VMN" = logFC) %>% unique(.) %>% arrange(-`F-VMN`)
m.vmn.g <- m_vmn %>% filter(genes %in% g)%>% select(genes, "M-VMN" = logFC) %>% unique(.) %>% arrange(-`M-VMN`)

fe.pit.g <- fe_pit %>% filter(genes %in% g)%>% select(genes, "F-PIT" = logFC) %>% unique(.) %>% arrange(-`F-PIT`)
m.pit.g <- m_pit %>% filter(genes %in% g)%>% select(genes, "M-PIT" = logFC) %>% unique(.) %>% arrange(-`M-PIT`)

ovt.g <- fe_ovt %>% filter(genes %in% g)%>% select(genes, "Ovaries" = logFC) %>% unique(.) %>% arrange(-Ovaries)
tes.g <- m_tes %>% filter(genes %in% g)%>% select(genes, "Testes" = logFC) %>% unique(.) %>% arrange(-`Testes`)

hm_g <-m.arc.g %>% full_join(fe.arc.g)  %>% full_join(m.vmn.g)%>% full_join(fe.vmn.g) %>% full_join(m.pit.g) %>% full_join(fe.pit.g)  %>% full_join(tes.g)%>% full_join(ovt.g) 
colnames(hm_g)
hm_gx <- hm_g %>% pivot_longer(cols = 2:9, names_to = "tissue", values_to = "NXM vs. Control")
# hm_gx[is.na(
library(ComplexHeatmap)
hm_gxx <- hm_g%>% arrange(genes) %>% column_to_rownames(., var= "genes") %>% as.matrix()
hm_gxx[is.na(hm_gxx)] <- 0

Heatmap(hm_gxx, name = "LogFC",cluster_columns = FALSE,cluster_rows = FALSE,
        show_heatmap_legend = T, column_names_rot = 60, rect_gp = gpar(col = "white", lwd = 2))







hm_all <- hm_c %>% rbind(hm_esr, hm_g, hm_im, hm_m, hm_s)

library(ComplexHeatmap)
hm_allxx <- hm_all%>% column_to_rownames(., var= "genes") %>% as.matrix()
hm_allxx[is.na(hm_allxx)] <- 0

hm_allxxT <- t(hm_allxx)
Heatmap(hm_allxxT, name = "LogFC",cluster_columns = FALSE,cluster_rows = F,
        show_heatmap_legend = T, column_names_rot = 60, rect_gp = gpar(col = "white", lwd = 2))

png("imgs/Heatmap_redo.png",width=17,height=3,units="in",res=1200)

Heatmap(hm_allxxT, name = "LogFC",cluster_columns = FALSE,cluster_rows = F,
        show_heatmap_legend = T, column_names_rot = 75, rect_gp = gpar(col = "white", lwd = 2))



dev.off()


hm_allx <- hm_c %>% rbind(hm_esr, hm_g)
library(ComplexHeatmap)
hm_allxx <- hm_allx%>% column_to_rownames(., var= "genes") %>% as.matrix()
hm_allxx[is.na(hm_allxx)] <- 0

hm_allxxT <- t(hm_allxx)
Heatmap(hm_allxxT, name = "LogFC",cluster_columns = FALSE,cluster_rows = F,
        show_heatmap_legend = T, column_names_rot = 45, rect_gp = gpar(col = "white", lwd = 2))

png("imgs/Heatmap_redo_partA.png",width=8.5,height=2.75,units="in",res=1200)

Heatmap(hm_allxxT, name = "LogFC",cluster_columns = FALSE,cluster_rows = F,
        show_heatmap_legend = T, column_names_rot =77, rect_gp = gpar(col = "white", lwd = 2))


dev.off()



hm_ally <-  hm_im %>% rbind(hm_m, hm_s)
library(ComplexHeatmap)
hm_allyy <- hm_ally%>% column_to_rownames(., var= "genes") %>% as.matrix()
hm_allyy[is.na(hm_allyy)] <- 0

hm_allyyT <- t(hm_allyy)
Heatmap(hm_allyyT, name = "LogFC",cluster_columns = FALSE,cluster_rows = F,
        show_heatmap_legend = T, column_names_rot = 45, rect_gp = gpar(col = "white", lwd = 2))

png("imgs/Heatmap_redo_partB.png",width=8.5,height=2.75,units="in",res=1200)

Heatmap(hm_allyyT, name = "LogFC",cluster_columns = FALSE,cluster_rows = F,
        show_heatmap_legend = T, column_names_rot =77, rect_gp = gpar(col = "white", lwd = 2))


dev.off()

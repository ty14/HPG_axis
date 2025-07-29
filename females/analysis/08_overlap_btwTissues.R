library(tidyverse)
my_logFC_threshold = 0.2

# females vmn
fe_vmn <- readRDS("females/clean_data/limma_eFDR_VMN.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val)  %>% distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold) %>%
  filter(.,P.Value <0.05)

vf_up <- fe_vmn %>% filter(logFC >= 0.2) %>% arrange(-logFC) # 138
vf_down <- fe_vmn %>% filter(logFC <= -0.2)%>% arrange(logFC) #131

# females arc
fe_arc <- readRDS("females/clean_data/limma_eFDR_ARC.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val) %>% distinct(.)  %>% 
  filter(.,abs(logFC) >= my_logFC_threshold) %>%
  filter(.,P.Value <0.05)
af_up <- fe_arc %>% filter(logFC >= 0.2) %>% arrange(-logFC) #79
af_down<- fe_arc %>% filter(logFC <= -0.2)%>% arrange(logFC) #113

# females pit
fe_pit <- readRDS("females/clean_data/limma_eFDR_PIT.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val)%>% distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold) %>%
  filter(.,P.Value <0.05)

pf_up <- fe_pit %>% filter(logFC >= 0.2) %>% arrange(-logFC) #110
pf_down <- fe_pit %>% filter(logFC <= -0.2)%>% arrange(logFC) #101


# females ovt
fe_ovt <- readRDS("females/clean_data/limma_eFDR_OVT.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val) %>% distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold) %>%
  filter(.,P.Value <0.05)

o_up <- fe_ovt %>% filter(logFC >= 0.2) %>% arrange(-logFC) #762
o_down <- fe_ovt %>% filter(logFC <= -0.2)%>% arrange(logFC) #533



#number of gene cut offs:
### Number of genes in each cut off 

#VMN
fe_vmn%>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #138
fe_vmn %>% filter(between(logFC, .2, .5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 121
fe_vmn %>% filter(between(logFC, .5, 1))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #15
fe_vmn %>% filter(between(logFC, 1, 8))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #2

fe_vmn %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 131 
fe_vmn %>% filter(between(logFC, -0.5, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #102
fe_vmn %>% filter(between(logFC, -1, -0.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 12
fe_vmn %>% filter(between(logFC, -8, -1)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 17


#ARC
fe_arc%>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #79
fe_arc %>% filter(between(logFC, .2, .5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 75
fe_arc %>% filter(between(logFC, .5, 1))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #4 
fe_arc %>% filter(between(logFC, 1, 8))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #0

fe_arc %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 113 
fe_arc %>% filter(between(logFC, -0.5, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #109
fe_arc %>% filter(between(logFC, -1, -0.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 3
fe_arc %>% filter(between(logFC, -8, -1)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 1


#PIT
fe_pit%>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #110
fe_pit %>% filter(between(logFC, .2, .5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 104
fe_pit %>% filter(between(logFC, .5, 1))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #6
fe_pit %>% filter(between(logFC, 1, 8))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #0

fe_pit %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 101
fe_pit %>% filter(between(logFC, -0.5, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #96
fe_pit %>% filter(between(logFC, -1, -0.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 4
fe_pit %>% filter(between(logFC, -8, -1)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 1



#OVT
fe_ovt%>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #762
fe_ovt %>% filter(between(logFC, .2, .5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 588
fe_ovt %>% filter(between(logFC, .5, 1))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #160
fe_ovt %>% filter(between(logFC, 1, 8))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #14

fe_ovt %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 533
fe_ovt %>% filter(between(logFC, -0.5, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #482
fe_ovt %>% filter(between(logFC, -1, -0.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 47
fe_ovt %>% filter(between(logFC, -8, -1)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 4




# upsetter Plot 
library(UpSetR)
library(workflowr)
library(ComplexUpset)
#VMN
listInput <- list(o_up= o_up$genes, v_up=v_up$genes,
                  o_down = o_down$genes, v_down = v_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

o_up$genes[o_up$genes %in% v_up$genes] 
#[1] "Kmt5a"    "Zfp709l1" "Wnk4"     "Susd2"    "Irgq"     "Fam107a"
o_down$genes[o_down$genes %in% v_down$genes]
# [1] "Col11a1" "Coq8a"  
o_down$genes[o_down$genes %in% v_up$genes] 
# [1] "Setbp1"   "Serpinb9" "Rhoj"     "Cd302"    "Cebpd"    "Slc16a9"  "Cyba"     "Htra3"    "Olfml3"  
o_up$genes[o_up$genes %in% v_down$genes]
# [1] "Ccdc78"     "Ccdc180"    "Lrrc34"     "Spef1"      "Cfap141"    "Rskr"       "LOC502684"  "Slc25a14"  
# [9] "RGD1565685" "Ccdc39"     "Cfap161"    "Cimap1b"    "Dnah10"     "T2"         "Krt8"       "Cfap251"   
# [17] "Ccdc153"    "Lrrc51"     "Rimklb"     "Hook2"      "Med15"      "Traf3ip1"   "Drc3"    

# ARC
listInput <- list(o_up= o_up$genes, a_up=a_up$genes,
                  o_down = o_down$genes, a_down = a_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

o_up$genes[o_up$genes %in% a_up$genes] 
# [1] "Itgb4"   "Ankzf1"  "Fam107a" "Snapc4"  "Dut" 
o_down$genes[o_down$genes %in% a_down$genes] 
# [1] "Nudt22"
o_down$genes[o_down$genes %in% a_up$genes]
# [1] "Sult5a1" "Dkk3" 
o_up$genes[o_up$genes %in% a_down$genes]
# [1] "Ppip5k1" "Gga3"    "Snrpd3"  "Trit1"  

#Pit
listInput <- list(o_up= o_up$genes, p_up=p_up$genes,
                  o_down = o_down$genes, p_down = p_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

o_up$genes[o_up$genes %in% p_up$genes] 
# [1] "Kmt5a"   "Meltf"   "Celsr2"  "Mcu"     "Irgq"    "mt-Rnr2"
o_down$genes[o_down$genes %in% p_down$genes] 
# [1] "Josd1"      "RT1-Db1"    "Dpt"        "Cd74"       "Snta1"      "Icoslg"     "Mix23"      "Emp3"      
# [9] "RGD1311164" "Lmo2"       "Cox4i2"     "Bicc1"
o_down$genes[o_down$genes %in% p_up$genes]
# [1] "Ntrk2"     "Lrch2"     "Serpinb9"  "LOC683456" "Csrp2"  
o_up$genes[o_up$genes %in% p_down$genes]
# [1] "Lmbr1l"       "P2ry4"        "LOC120095166" "Cenatac" "Syngap1"      "Ints6l"       "Zscan12"    



mat1 <- matrix(c(6,9,23,2),nrow=2)
chisq.test(mat1)
psych::phi(mat1)

mat2 <- matrix(c(5,6,4,1),nrow=2)
chisq.test(mat2)
psych::phi(mat2)


mat3 <- matrix(c(6,5,7,12),nrow=2)
chisq.test(mat3)
psych::phi(mat3)



listInput <- list(a_up= a_up$genes, v_up=v_up$genes,
                  a_down = a_down$genes, v_down = v_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

v_up$genes[v_up$genes %in% a_up$genes] 
a_down$genes[a_down$genes %in% v_down$genes] 
v_up$genes[v_up$genes %in% a_down$genes]
a_up$genes[a_up$genes %in% v_down$genes]


listInput <- list(a_up= a_up$genes, p_up=p_up$genes,
                  a_down = a_down$genes, p_down = p_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

p_up$genes[p_up$genes %in% a_up$genes] 
a_down$genes[a_down$genes %in% p_down$genes] 
p_up$genes[p_up$genes %in% a_down$genes]


listInput <- list(v_up= v_up$genes, p_up=p_up$genes,
                  v_down = v_down$genes, p_down = p_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

v_up$genes[v_up$genes %in% p_up$genes] 

p_up$genes[p_up$genes %in% v_down$genes]

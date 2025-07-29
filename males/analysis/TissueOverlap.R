library(tidyverse)
my_logFC_threshold = 0.2

# males vmn
vmn <- readRDS("males/clean_data/limma_eFDR_VMN.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val)  %>% distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold) %>%
  filter(.,P.Value <0.05)

v_up <- vmn %>% filter(logFC >= 0.2) %>% arrange(-logFC) # 118
v_down <- vmn %>% filter(logFC <= -0.2)%>% arrange(logFC) #85

# males arc
arc <- readRDS("males/clean_data/limma_eFDR_ARC.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val) %>% distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold) %>%
  filter(.,P.Value <0.05)
a_up <- arc %>% filter(logFC >= 0.2) %>% arrange(-logFC) #317
a_down<- arc %>% filter(logFC <= -0.2)%>% arrange(logFC) #296

# males pit
pit <- readRDS("males/clean_data/limma_eFDR_PIT.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val)%>% distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold) %>%
  filter(.,P.Value <0.05)

p_up <- pit %>% filter(logFC >= 0.2) %>% arrange(-logFC) #84
p_down <- pit %>% filter(logFC <= -0.2)%>% arrange(logFC) #110


# males test
tes <- readRDS("males/clean_data/limma_eFDR_Testis.RDS") %>% 
  dplyr::select(genes, logFC, P.Value,  adj.P.Val) %>% distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold) %>%
  filter(.,P.Value <0.05)

t_up <-tes%>% filter(logFC >= 0.2) %>% arrange(-logFC) #65
t_down <-tes %>% filter(logFC <= -0.2)%>% arrange(logFC) #88



#number of gene cut offs:
### Number of genes in each cut off 

#VMN
vmn%>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #118
vmn %>% filter(between(logFC, .2, .5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 102
vmn %>% filter(between(logFC, .5, 1))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #6
vmn %>% filter(between(logFC, 1, 8))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #10

vmn %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 85 
vmn %>% filter(between(logFC, -0.5, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #75
vmn %>% filter(between(logFC, -1, -0.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 6
vmn %>% filter(between(logFC, -8, -1)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 4


#ARC
arc%>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #317
arc %>% filter(between(logFC, .2, .5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 241
arc %>% filter(between(logFC, .5, 1))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #69
arc %>% filter(between(logFC, 1, 8))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #7

arc %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 296
arc %>% filter(between(logFC, -0.5, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #172
arc %>% filter(between(logFC, -1, -0.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 88
arc %>% filter(between(logFC, -8, -1)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 36


#PIT
pit%>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #84
pit %>% filter(between(logFC, .2, .5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 77
pit %>% filter(between(logFC, .5, 1))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #6
pit %>% filter(between(logFC, 1, 8))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #1

pit %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 110
pit %>% filter(between(logFC, -0.5, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #106
pit %>% filter(between(logFC, -1, -0.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 4
pit %>% filter(between(logFC, -8, -1)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 0



#tes
tes%>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #152
tes %>% filter(between(logFC, .2, .5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 149
tes %>% filter(between(logFC, .5, 1))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #2
tes %>% filter(between(logFC, 1, 8))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #1

tes %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 116
tes %>% filter(between(logFC, -0.5, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #108
tes %>% filter(between(logFC, -1, -0.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 7
tes %>% filter(between(logFC, -8, -1)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 1




# upsetter Plot 
library(UpSetR)
library(workflowr)
library(ComplexUpset)
#VMN
listInput <- list(t_up= t_up$genes, v_up=v_up$genes,
                  t_down = t_down$genes, v_down = v_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

t_up$genes[t_up$genes %in% v_up$genes] 
#[1] "Sphkap"
t_down$genes[t_down$genes %in% v_down$genes]
# [1] "Gas6" 

# ARC
listInput <- list(t_up= t_up$genes, a_up=a_up$genes,
                  t_down = t_down$genes, a_down = a_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

t_up$genes[t_up$genes %in% a_up$genes] 
# [1] "Stk16" 
t_down$genes[t_down$genes %in% a_down$genes] 
# [1] "Col1a1" "Igfbp6" "Ogn"
t_up$genes[t_up$genes %in% a_down$genes]
# [[1] "Sfi1" 

#Pit
listInput <- list(t_up= t_up$genes, p_up=p_up$genes,
                  t_down = t_down$genes, p_down = p_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

t_up$genes[t_up$genes %in% p_up$genes] 
# [1] "Inip"
t_down$genes[t_down$genes %in% p_down$genes] 
# [1] "Gsdma"


mat1 <- matrix(c(21,0,0,6),nrow=2)
chisq.test(mat1)
psych::phi(mat1)

mat2 <- matrix(c(5,0,1,1),nrow=2)
chisq.test(mat2)
psych::phi(mat2)


mat3 <- matrix(c(2,1,0,3),nrow=2)
chisq.test(mat3)
psych::phi(mat3)



listInput <- list(a_up= a_up$genes, v_up=v_up$genes,
                  a_down = a_down$genes, v_down = v_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

v_up$genes[v_up$genes %in% a_up$genes] 
a_down$genes[a_down$genes %in% v_down$genes] 
# > v_up$genes[v_up$genes %in% a_up$genes] 
# [1] "Rprm"         "Gpr83"        "Akap5"        "Me3"          "Prkcg"        "Nefm"        
# [7] "Kcnj11"       "Reps2"        "Cdk14"        "Add2"         "Lingo1"       "Neto1"       
# [13] "Mgst3"        "Cox10"        "Ccdc184"      "Bdnf"         "Dgkb"         "Cacnb4"      
# [19] "Tmem14a"      "Ppfia2"       "LOC108353239"
# > a_down$genes[a_down$genes %in% v_down$genes] 
# [1] "Gstm2"   "Flnc"    "Vamp8"   "Sdc2"    "Heatr6"  "Dennd2b"


listInput <- list(a_up= a_up$genes, p_up=p_up$genes,
                  a_down = a_down$genes, p_down = p_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

p_up$genes[p_up$genes %in% a_up$genes] 
a_down$genes[a_down$genes %in% p_down$genes] 
a_up$genes[a_up$genes %in% p_down$genes]
# [1] "Cacnb4" "Riox1"  "Lin7a"  "Pdp1"   "Zfp248"
# > a_down$genes[a_down$genes %in% p_down$genes] 
# [1] "Txnip"
# > a_up$genes[a_up$genes %in% p_down$genes]
# [1] "Trpc5"

listInput <- list(v_up= v_up$genes, p_up=p_up$genes,
                  v_down = v_down$genes, p_down = p_down$genes)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

v_up$genes[v_up$genes %in% p_up$genes] 
p_down$genes[p_down$genes %in% v_down$genes]
p_up$genes[p_up$genes %in% v_down$genes]
> v_up$genes[v_up$genes %in% p_up$genes] 
# [1] "Cacnb4" "Adora1"
# > p_down$genes[p_down$genes %in% v_down$genes]
# [1] "Slc7a6" "Apex2"  "Rpp40" 
# > p_up$genes[p_up$genes %in% v_down$genes]
# [1] "Hacl1"
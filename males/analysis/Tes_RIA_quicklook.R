library(tidyverse)

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



p1 <- ggplot(nmx, aes(x = treatment, y = con_ng_ml)) + 
  geom_boxplot(aes( fill = treatment), outlier.shape = NA) +
  geom_jitter(width=0.15, color="black", alpha = 1) +
  theme_classic() +
  ggtitle("with high CV sample")

p1

p2 <- ggplot(m_data, aes(x = treatment, y = con_ng_ml)) + 
  geom_boxplot(aes( fill = treatment), outlier.shape = NA) +
  geom_jitter(width=0.15, color="black", alpha = 1) +
  scale_color_manual(values = c( "darkorange", "blue"))+
  scale_fill_manual(values = c( "darkorange", "blue"))+    
  theme_classic() +
  theme(legend.position = "none")

p2

gridExtra::grid.arrange(p1, p2, ncol = 2)


library(lmerTest)
hist(m_data$con_ng_ml)

m1 <- lmer(con_ng_ml ~ treatment+ (1|Cohort), data = m_data)
summary(m1)


m1 <- glmer(round(con_ng_ml) ~ treatment+ (1|Cohort),family=poisson, data = m_data)
summary(m1)

sf <- df2 %>% filter(!SampleName %in% nmx$SampleName)
sf <- sf[c(11:76),]
sf$SampleName


p3 <- ggplot(sf, aes(x = treatment, y = con_ng_ml)) + 
  geom_boxplot(aes( fill = treatment), outlier.shape = NA) +
  geom_jitter(width=0.15, color="black", alpha = 1) +
  theme_classic() +
  ggtitle("with high CV sample")

p3

sf2 <- sf %>% filter(cv <20)

p4 <- ggplot(sf2, aes(x = treatment, y = con_ng_ml)) + 
  geom_boxplot(aes( fill = treatment), outlier.shape = NA) +
  geom_jitter(width=0.15, color="black", alpha = 1) +
  theme_classic() +
  ggtitle("without high CV samples")

p4

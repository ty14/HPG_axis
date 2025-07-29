# Assuming "expression_data" is your gene expression matrix

library(limma)
library(edgeR)
library(tidyverse)


# sample information 
id <- read_csv("raw_data/traitData.csv")
colnames(id)
idx <- id %>% select(1:4,15)

##### ARC
#counts 
a <- readRDS("brain/clean_data/arc_counts.RDS") %>%
  select(-region) %>% 
  column_to_rownames(., "ensgene")

colnames(a) <- paste("A", colnames(a), sep = "_")

a_data <- idx
a_data$SampleName <- paste("A", a_data$SampleName, sep = "_")
a_datax <- a_data %>% column_to_rownames(., var = "SampleName")


#checks
all(row.names(a_datax) %in% colnames(a)) #check 
a_datax <- a_datax[colnames(a),]
a <- a[,rownames(a_datax)]
all(rownames(a_datax) == colnames(a)) #check

#just females
fa_datax <- a_datax %>% filter(Sex == "F")
all(row.names(fa_datax) %in% colnames(a)) #check 
a <- a[,rownames(fa_datax)]
all(rownames(fa_datax) == colnames(a)) #check

mostVar <- function(data, n, i_want_most_var = TRUE) {
  data.var <- apply(data, 1, stats::var)
  data[order(data.var, decreasing = i_want_most_var)[1:n],] 
}


arc_var <- mostVar(a, 4000, i_want_most_var = TRUE)

##### VMN
#counts 
v <- readRDS("brain/clean_data/VMN_counts.RDS") %>%
  select(-region) %>% 
  column_to_rownames(., "ensgene")

colnames(v) <- paste("V", colnames(v), sep = "_")

v_data <- idx
v_data$SampleName <- paste("V", v_data$SampleName, sep = "_")
v_datax <- v_data %>% column_to_rownames(., var = "SampleName")


#checks
all(row.names(v_datax) %in% colnames(v)) #check 
v_datax <- v_datax[colnames(v),]
v <- v[,rownames(v_datax)]
all(rownames(v_datax) == colnames(v)) #check



#just females
fv_datax <- v_datax %>% filter(Sex == "F")
all(row.names(fv_datax) %in% colnames(v)) #check 
v <- v[,rownames(fv_datax)]
all(rownames(fv_datax) == colnames(v)) #check

mostVar <- function(data, n, i_want_most_var = TRUE) {
  data.var <- apply(data, 1, stats::var)
  data[order(data.var, decreasing = i_want_most_var)[1:n],] 
}


vmn_var <- mostVar(v, 4000, i_want_most_var = TRUE)


##### PIT
#counts 
p <- readRDS("brain/clean_data/PIT_counts.RDS") %>%
  select(-region) %>% 
  column_to_rownames(., "ensgene")

colnames(p) <- paste("P", colnames(p), sep = "_")

p_data <- idx
p_data$SampleName <- paste("P", p_data$SampleName, sep = "_")
p_datax <- p_data %>% column_to_rownames(., var = "SampleName")


#checks
all(row.names(p_datax) %in% colnames(p)) #check 
p_datax <- p_datax[colnames(p),]
p <- p[,rownames(p_datax)]
all(rownames(p_datax) == colnames(p)) #check

#just females
fp_datax <- p_datax %>% filter(Sex == "F")
all(row.names(fp_datax) %in% colnames(p)) #check 
p <- p[,rownames(fp_datax)]
all(rownames(fp_datax) == colnames(p)) #check

mostVar <- function(data, n, i_want_most_var = TRUE) {
  data.var <- apply(data, 1, stats::var)
  data[order(data.var, decreasing = i_want_most_var)[1:n],] 
}


pit_var <- mostVar(p, 4000, i_want_most_var = TRUE)


#OVT
df <- read.table("raw_data/all_gonad_counts_Star.gff", col.names = c("symbol", "count","id"))
head(df)
dfx <- df %>% pivot_wider(names_from = id, values_from = count) %>% column_to_rownames(., var = "symbol")
colnames(dfx) <- substr(colnames(dfx), 7, 12)
# dfxx <- dfx %>% rownames_to_column(., var = "symbol")
# saveRDS(dfxx, "raw_data/gonads_clean_counts.RDS")

#counts
g <- dfx


#coldata vmn
g_data <- idx
g_datax <- g_data %>% column_to_rownames(., var = "SampleName")
fg_datax <- g_datax %>% filter(Sex == "F")

#checks
all(row.names(fg_datax) %in% colnames(g)) #check 
g <- g[,rownames(fg_datax)]
all(rownames(fg_datax) == colnames(g)) #check


mostVar <- function(data, n, i_want_most_var = TRUE) {
  data.var <- apply(data, 1, stats::var)
  data[order(data.var, decreasing = i_want_most_var)[1:n],] 
}


ovt_var <- mostVar(g, 4000, i_want_most_var = TRUE)


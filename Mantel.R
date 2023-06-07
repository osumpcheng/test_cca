library("vegan")

df <- read.csv("Creek_order_0310.csv", header=TRUE, row.names = 1)
df_log <- log10(df+1)
df_log <- df_log[colSums(df_log)>0]
dist.abund = vegdist(df_log, method = "bray")
man_al <- rep(0,1000)

for (i in seq(1,1000)){
  df_50 <- sample(df_log, size = 50, replace = FALSE)
  dist_50 <- vegdist(df_50, method = "bray")
  mant <- mantel(dist.abund, dist_50, method = "spearman", permutations = 999, na.rm = TRUE)
  man_al[i] <- mant$statistic
  if (mant$statistic>=max(man_al)){
    df_51 <- df_50}
}

df_chem <- read.csv("Chem_0924_46.csv", header=TRUE, row.names = 1)
dist.chem = vegdist(df_chem, method = "bray")
man_chem_bio <- rep(0,1000)

for (i in seq(1,1000)){
  df_50 <- sample(df_log, size = 50, replace = FALSE)
  dist_50 <- vegdist(df_50, method = "bray")
  mant <- mantel(dist.chem, dist_50, method = "spearman", permutations = 999, na.rm = TRUE)
  man_chem_bio[i] <- mant$statistic
  if (mant$statistic>=max(man_chem_bio)){
    df_52 <- df_50}
}

mantel_0 <- mantel(dist.chem, dist.abund, method = "spearman", permutations = 999, na.rm = TRUE)

library("vegan")
library("ggplot2")
micro <- read.csv('creek_microb_otu.csv',row.names=1)
ord <- read.csv('creek_microb_order.csv',row.names=1)

# NMS by all w/o log transform
micro_mds <- metaMDS(micro)
plot(micro_mds,type ="n")
points(micro_mds, display = "sites", cex = 1, pch=0, col="red")
text(micro_mds, display = "sites", cex=1, col="blue")

micro_log <- log10(micro+1)
micro_logmds <- metaMDS(micro_log)
plot(micro_logmds,type ="n")
points(micro_logmds, display = "sites", cex = 1, pch=0, col="red")
text(micro_logmds, display = "sites", cex=0.5, col="blue")

# compare ordination
pro <- procrustes(micro_mds, micro_logmds)
plot(pro)

# NMS by order
ord_mds <- metaMDS(ord)
plot(ord_mds$points,type ="n")
points(ord_mds, display = "sites", cex = 1, pch=0, col="red")
text(ord_mds, display = "sites", cex=0.7, col="blue")

pro <- procrustes(micro_mds, ord_mds)
plot(pro)

# fit NMS by ord, find determinant ord, too crowded, no subset worked out
ord.fit <- envfit(ord_mds, ord, perm=999)
ord.fit
plot(ord.fit, col="gray50", cex=0.8, p.max=0.001, font=4)
#ord.fit$vectors$arrows[ord.fit$vectors$r>0.5,]

# plot groupings (ordihull)
env <- read.csv('sample_descrip.txt',row.names = 1, sep = "\t",check.names = FALSE)
env <- na.omit(t(env))
ordihull(micro_mds, env[,'Location'])

# remove WWTP data amd NMS with creek only (both all and order)
wwtp <- c("18B1A", "18B2A", "18B3A", "18B4A","18E1A", "18E2A", "18E3A", "18E4A", "18I1A", "18I2A", "18I3A", "18I4A", "19E1A", "19E2A", "19E3A", "19E4A")
env_c <- env[!(row.names(env) %in% wwtp), ]
micro_c <- micro[!(row.names(micro) %in% wwtp), ]
micro_c <- micro_c[, colSums(micro_c != 0) > 0]
c_mds <- metaMDS(micro_c)
plot(c_mds,type ="n")
points(c_mds, display = "sites", cex = 1, pch=0, col="black")
#text(c_mds, display = "sites", cex = 0.8, col="blue")
ordihull(c_mds, groups=env_c[,'Location'],label=TRUE, col=c('brown','forestgreen','purple','red'), lwd=2)


scor_site <- as.data.frame(scores(c_mds)$sites)
#scor_site$grp <- paste0(env_c[,2], env_c[,1])
scor_site$grp <- env_c[,2]
#scor_site$site <- rownames(chem_55)
#scor_site$crk <- as.factor()
#scor_site$sec <- as.character(env_55$Section)
 grp_spring <- scor_site[scor_site$grp == "Spring", ][chull(scor_site[scor_site$grp == 
                                                           "Spring", c("NMDS1", "NMDS2")]), ] 
 grp_fall <- scor_site[scor_site$grp == "Fall", ][chull(scor_site[scor_site$grp == 
                                                                        "Fall", c("NMDS1", "NMDS2")]), ] 
 grp_winter <- scor_site[scor_site$grp == "Winter", ][chull(scor_site[scor_site$grp == 
                                                                        "Winter", c("NMDS1", "NMDS2")]), ] 
 hull.data <- rbind(grp_spring, grp_fall, grp_winter)
ggplot()+
  #geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=grp,group=grp),alpha=0.30) + # add the convex hulls
  geom_point(data=scor_site, aes(x= NMDS1, y=NMDS2, shape=grp, fill= grp, colour=grp), size=5)+
  #geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=grp,group=grp),alpha=0.30) + # add the convex hulls
  #geom_text(data=scor_site,aes(x=NMDS1,y=NMDS2,label=site),size=4,vjust=0,hjust=0)+
  scale_shape_manual(values=c( 18, 16, 15, 17, 18, 16, 15, 17, 18, 16, 15))+
  scale_colour_manual(values=c('orange','orange','orange','green4','green4','green4','green4',
                               'deeppink', 'deeppink', 'deeppink','deeppink'))+
  theme_bw() +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())



ord_c <- ord[!(row.names(ord) %in% wwtp), ]
ord_c <- ord_c[, colSums(ord_c != 0) > 0]
o_mds <- metaMDS(ord_c)
plot(o_mds$points,type ="n")
points(o_mds, display = "sites", cex = 1, pch=0, col="red")
text(o_mds, display = "sites", cex = 0.8, col="blue")
ordihull(o_mds, env_c[,'Type'], col='black')

pro <- procrustes(c_mds, o_mds)
plot(pro)

c.fit <- envfit(c_mds, ord_c, perm=999)
c.fit
plot(c.fit, col="gray50", p.max=0.05, cex=0.8, font=4)


#Chemical NMDS
chem_55 <- read.csv('Chem_0924_55.csv',row.names = 1, numerals ="no.loss")
env_55 <- read.csv('EnvVar_nldas.csv',row.names = 1)

chem.mds <- metaMDS(chem_55, distance = "bray", k = 3, autotransform = FALSE, noshare = TRUE, trace = TRUE, plot = TRUE, wascores = TRUE, parallel =4)
chem_sco <- as.data.frame(chem.mds$points[,1:2])
chem_fit <- envfit(chem_sco, env_55[,4:25], perm=999)
chem_fit.df <- as.data.frame(chem_fit$vectors$arrows*sqrt(chem_fit$vectors$r))
chem_fit.df$species<-rownames(chem_fit.df)

ggplot(data = chem_sco, aes(x = MDS1, y = MDS2, label = rownames(chem_sco)))+
  geom_line(data = chem_sco, aes(group=paste0(env_55$Creek, env_55$Month)),  color="black", size=1, 
            alpha=0.5, arrow = arrow(length=unit(0.20,"cm"), ends="last", type = "closed")) +
  geom_point(data = chem_sco, aes(colour = env_55$Month), size = 4, alpha = 0.5) +
  scale_colour_manual(values = c("orange", "steelblue", "coral", "darkseagreen")) +
  geom_segment(data=chem_fit.df,aes(x=0,xend=MDS1,y=0,yend=MDS2), 
               arrow = arrow(length = unit(0.5, "cm")),colour="red",inherit_aes=FALSE) + 
  geom_text(data=chem_fit.df,aes(x=MDS1+0.01,y=MDS2+0.01,label=species),size=5)

  
## batch effect
batch_chem <- c(rep('batch_1', 20), rep('batch_2', 10), rep('batch_3', 10), rep('batch_4', 15))
library(sva)
#library(bladderbatch)
#data(bladderdata)
#dat_blad <- bladderEset[1:50,]
#pheno_blad = pData(dat_blad)
#edata = exprs(dat_blad)
combat_chem <- ComBat(dat=t(chem_55), batch=batch_chem, mod=NULL, par.prior=TRUE)
combat_chem.mds <- metaMDS(t(combat_chem), distance = "bray", k = 3, autotransform = FALSE, noshare = TRUE, trace = TRUE, plot = TRUE, wascores = TRUE, parallel =4)

library(VennDiagram)

env_55$Season <- c(rep('Summer', 20), rep('Winter', 20), rep('Fall', 15))
batch_all <- matrix(nrow = ncol(chem_55), ncol = 3)

for (i in c(1:3)){
  bat <- unique(env_55$Season)[i]
  bat_id <- which(env_55$Season==bat)
  print(row.names(chem_55[bat_id,]))
  batch_all[,i] <- apply(chem_55[bat_id,], 2, max)>0
}
venn_batch <- list(
  Summer = colnames(chem_55)[batch_all[,1]], 
  Winter = colnames(chem_55)[batch_all[,2]], 
  Fall = colnames(chem_55)[batch_all[,3]]
  #section_4 = colnames(chem_55)[batch_all[,4]],
  #section_5 = colnames(chem_55)[batch_all[,5]]
)
venn.diagram(venn_batch, filename = "venn-4-seasons.png")
chem_pre <- chem_55
chem_pre[chem_pre>0] <- 1
over1_chem <- chem_55[,rowSums(t(chem_pre))>1]
over4_chem <- chem_55[,rowSums(t(chem_pre))>4]
for (i in c(1:4)){
  bat <- unique(batch_chem)[i]
  assign(bat, colnames(over1_chem)[apply(over1_chem[which(batch_chem==bat),], 2, max)>0])
}
venn_batch <- list(batch_1=batch_1, batch_2=batch_2, batch_3=batch_3, batch_4=batch_4)
venn.diagram(venn_batch, filename = "venn-4-batches_over1.png")

## install plsda batch
cran.pkgs <- c('pheatmap', 'vegan', 'ruv', 'UpSetR', 'gplots', 
               'ggplot2', 'gridExtra', 'performance', 'BiocManager')

for(c in seq_len(length(cran.pkgs))){
  if (!requireNamespace(cran.pkgs[c], quietly = TRUE))
    install.packages(cran.pkgs[c])
}

# Bioconductor
bioc.pkgs <- c('mixOmics', 'sva', 'limma', 'Biobase', 'metagenomeSeq')

for(b in seq_len(length(bioc.pkgs))){
  if (!requireNamespace(bioc.pkgs[b], quietly = TRUE))
    BiocManager::install(bioc.pkgs[b])
}
library(PLSDAbatch)
#data('AD_data')
combat_chem_pl <- PLSDA_batch(X=chem_55, Y.bat=batch_chem, ncomp.trt = 1, ncomp.bat = 4, balance = FALSE)
plbat_chem.mds <- metaMDS(combat_chem_pl$X.nobatch, distance = "bray", k = 3, autotransform = FALSE, noshare = TRUE, trace = TRUE, plot = TRUE, wascores = TRUE, parallel =4)
over4_mds <- metaMDS(over4_chem, distance = "bray", k = 3, autotransform = FALSE, noshare = TRUE, trace = TRUE, plot = TRUE, wascores = TRUE, parallel =4)

scor_site <- as.data.frame(chem.mds$points) #chem.poa$points
scor_site <- as.data.frame(scores(chem.mds)$sites)
scor_site$grp <- env_55$Month
scor_site$site <- rownames(chem_55)
scor_site$bat <- batch_chem
scor_site$crk <- env_55$Creek
ggplot()+
  geom_point(data=scor_site, aes(x= NMDS1, y=NMDS2, colour=bat), size=5)+
  geom_text(data=scor_site,aes(x=NMDS1,y=NMDS2,label=site),size=4,vjust=0,hjust=0)+
  scale_colour_manual(values=RColorBrewer::brewer.pal(4, "Set1"))


chem_norm <- read.csv('Normalized creeks chem.csv')
chem_mzrt <- chem_norm[,1]
chem_norm <- as.data.frame(chem_norm[,2:ncol(chem_norm)])
chem_norm <- chem_norm[apply(chem_norm, 1, sum)>1,]
chem_norm <- log10(chem_norm+1)

norm_chem.mds <- metaMDS(t(chem_norm), distance = "bray", k = 3, pc = TRUE, autotransform = TRUE,
                         noshare = TRUE, trace = TRUE, wascores = TRUE, parallel =4)


chem_all_raw <- read.csv('C:/Users/shiche/Documents/R/test_cca/normalized peaktable.txt', sep = "\t", skip=2)
chem_all_raw.loc <- chem_all_raw[1,]
colnames(chem_all_raw) <- colnames(chem_all_raw.loc)
batch1_raw <- chem_all_raw[,1:125]
batch2_raw <- chem_all_raw[,126:189]
batch3_raw <- chem_all_raw[,190:325]
batch4_raw <- chem_all_raw[,330:420]

batch1_grp <- c('18B', '18C', '18D', '18E', '18F', '18G')
batch2_grp <- c('19A', '19B')
batch3_grp <- c('19C', '19D', '19E')
batch4_grp <- c('18I', '18L', '18M', '18N')
wwtp <- c('18B', '18E', '18I', '19E')
batch4_raw <- batch4_raw[,c(5,c(14:82), 91)]
batch3_raw <- batch3_raw[,c(1, 9, c(11:46), c(67:76), c(80:91))]
batch2_raw <- batch2_raw[,c(1, 9, c(17:53), 63, 64)]
batch1_raw <- batch1_raw[,c(2, c(11:114), c(123:125))]


batch1_sample <- c(paste0(batch1_grp[1], c(1:4)), 
                   paste0(batch1_grp[2], c(1:5)), 
                   paste0(batch1_grp[3], c(1:5)),
                   paste0(batch1_grp[4], c(1:4)), 
                   paste0(batch1_grp[5], c(1:5)), 
                   paste0(batch1_grp[6], c(1:5)))
batch2_sample <- c(paste0(batch2_grp[1], c(1:5)), 
                   paste0(batch2_grp[2], c(1:5)))
batch3_sample <- c(paste0(batch3_grp[1], c(1:5)), 
                   paste0(batch3_grp[2], c(1:5)), 
                   paste0(batch3_grp[3], c(1:5)))
batch4_sample <- c(paste0(batch4_grp[1], c(1:4)), 
                   paste0(batch4_grp[2], c(1:5)), 
                   paste0(batch4_grp[3], c(1:5)),
                   paste0(batch4_grp[4], c(1:5)))


# merge triplicates
sample_name <- batch3_sample
batch_name <- batch3_grp
batch_raw <- batch3_raw
mergetable_3 <- matrix(nrow = nrow(batch_raw), ncol = length(sample_name))
colnames(mergetable_3) <- sample_name

for (grp_name in batch_name) {
  idx <- grep(grp_name, sample_name)
  raw_col <- grep(grp_name, colnames(batch_raw))
  ref_blk <- batch_raw[,c(raw_col[1], raw_col[1]-1, raw_col[1]-2, 
                       raw_col[length(raw_col)]+1, raw_col[length(raw_col)]+2)]
  #if (grp_name == '19A'){
  #  blank_col <- c(1,2,3,25,26)}
  #else{
  #  blank_col <- c(21,25,26,39,40)
  #}
  #ref_blk <- batch_raw[,blank_col]
  print(colnames(ref_blk))
  ref_all <- apply(ref_blk, 1, FUN = max)
  for (i in idx){
    grp_peaks <- batch_raw[,grep(sample_name[i], colnames(batch_raw))]
    print(colnames(grp_peaks))
    min_sample <- apply(grp_peaks, 1, FUN = min)
    print('min_sample')
    print(length(which(min_sample>100)))
    print('max_blank')
    print(length(which(ref_all>100)))
    real <- min_sample/ref_all >3 & min_sample > 100
    print(length(which(real==TRUE)))
    mergetable_4[real, i] <- rowMeans(grp_peaks[real,])
  } 
}

merged_raw = rbind(t(mergetable_1), t(mergetable_2), t(mergetable_3), t(mergetable_4))
colnames(merged_raw)<- chem_all_raw[,1]
merged.peaktable <- merged_raw[,colSums(is.na(merged_raw)) != nrow(merged_raw)]
merged.peaktable <- merged.peaktable[,-grep('group', colnames(merged.peaktable))]
chem.count <- rowSums(is.na(merged.peaktable)==FALSE)
merged.creek <- merged.peaktable[c(c(5:14), c(19:48), c(58:72)),]
merged.creek <- merged.creek[,colSums(is.na(merged.creek)) != nrow(merged.creek)]
merged.creek[is.na(merged.creek)] <- 0

count_55 <- rowSums(chem_55>0)

norm_mds <- metaMDS(merged.creek, distance = "bray", k = 3, autotransform = FALSE, 
                    noshare = TRUE, trace = TRUE, plot = TRUE, wascores = TRUE, parallel =4)

scor_site <- as.data.frame(scores(norm_mds)$sites)
scor_site$grp <- env_55$Month
scor_site$site <- rownames(chem_55)
scor_site$bat <- batch_chem
scor_site$crk <- env_55$Creek
scor_site$sec <- as.character(env_55$Section)
ggplot()+
  geom_point(data=scor_site, aes(x= NMDS1, y=NMDS2, colour=grp), size=5)+
  geom_text(data=scor_site,aes(x=NMDS1,y=NMDS2,label=site),size=4,vjust=0,hjust=0)+
  scale_colour_manual(values=RColorBrewer::brewer.pal(5, "Set1"))


cor_axs <- read.csv('chem_nms_3axes.csv', row.names = 1)
cor_env <- read.csv('EnvVar_nldas_55_renew.csv', row.names = 1)

for (i in c(1:22)){
  print(colnames(cor_env)[i])
  print(cor.test(cor_axs[,3], cor_env[,i]))
}

pca_sem_df <- read.csv('sem_significant_variables_46.csv', row.names = 1)
microb_axes <- read.csv('microb_nms_2axes.csv', row.names = 1)
chem_axes_55 <- read.csv('chemical_NMS.csv', row.names = 1)
chem_raw_log <- read.csv('C:/Users/shiche/Documents/R/test_cca/norm_merg_chem.csv', row.names = 1)
chem_raw_log <- log10(chem_raw_log+1)

cor.test(rowMeans(chem_raw_log), chem_axes_55[,2])
cor.test(rowSums(chem_raw_log!=0), chem_axes_55[,1])

trans_std <- function(data.week){
  return((data.week - mean(data.week, na.rm = TRUE))/sd(data.week, na.rm = TRUE))
}

pca_env <- prcomp(pca_sem_df[,1:8], scale. = TRUE)
pca_lc <- prcomp(pca_sem_df[,9:15])

trans_std <- function(data.week){
  return(data.week /sd(data.week, na.rm = TRUE))
}
pca_env_df <- pca_sem_df[,1:8]
pca_env_df[,5] <- trans_std(pca_env_df[,5])*10
pca_env_df[,4] <- trans_std(pca_env_df[,4])*10
pca_env_df[,3] <- pca_env_df[,3]/10
pca_env_df[,7] <- pca_env_df[,7]/10
pca_env_1 <- prcomp(pca_env_df[,1:8])

library(lavaan)
library(readxl)
library(semPlot)

sem_creeks <- cbind(pca_env$x[,1:3], pca_lc$x[,1:2], microb_axes, pca_sem_df[,16:18])
colnames(sem_creeks) <- c('Env_PC1', 'Env_PC2', 'Env_PC3', 'LC_PC1', 'LC_PC2', 'Micro_x1'
                          , 'Micro_x2', 'Chem_x1', 'Chem_x2', 'Chem_x3')

m1b_crk <- '
  #regressions
  Micro_x1 ~ Env_PC1 + Env_PC2 + Env_PC3 + LC_PC1
  Micro_x2 ~ Env_PC1 + Env_PC2
  Chem_x1 ~ Env_PC1 + Env_PC2 + Env_PC3 + Micro_x2
  Chem_x2 ~ Env_PC1 + Env_PC2 + LC_PC1 + LC_PC2 + Micro_x1
  Chem_x3 ~ Env_PC1 + Env_PC3 + LC_PC1 + LC_PC2 + Micro_x2
'
m2b_crk <- '
  #regressions
  Micro_x1 ~ a1*Env_PC1 + b1*Env_PC2 + c1*Env_PC3
  Micro_x2 ~ a2*Env_PC1 + f2*Chem_x1
  Chem_x1 ~ b3*Env_PC2 + c3*Env_PC3
  Chem_x2 ~ a4*Env_PC1 + b4*Env_PC2 + g4*LC_PC1 + h4*LC_PC2 + d4*Micro_x1
  Chem_x3 ~ b5*Env_PC2 + c5*Env_PC3 + e5*Micro_x2
  
  #mediated
  E2_C1_M2 := b3*f2
  E3_C1_M2 := c3*f2
  E1_M1_C2 := a1*d4
  E2_M1_C2 := b1*d4
  E3_M1_C2 := c1*d4
  E1_M2_C3 := a2*e5
  C1_M2_C3 := f2*e5
'

m3b_crk <- '
  #regressions
  Micro_x1 ~ Env_PC1 + Env_PC2 + Env_PC3
  Micro_x2 ~ Env_PC1 + Env_PC2
  Chem_x1 ~ Env_PC2 + Env_PC3 + Micro_x1 + Micro_x2
  Chem_x2 ~ Env_PC1 + Env_PC2 + LC_PC1 + LC_PC2 + Micro_x1
  Chem_x3 ~ Env_PC1 + Env_PC3 + LC_PC1 + LC_PC2 + Micro_x2
'

fit_crk <- sem(
  model          = m3b_crk, 
  data           = sem_creeks
)
#semPaths(fit_crk, "par", nCharNodes = 0, style = "lisrel", 
#         edge.label.cex = 1.2, rotation = 1,intercepts = FALSE, layout = 'circle2')
semPaths(fit_crk, whatLabels = "std", style = "ram", layout = "spring")
summary(fit_crk, fit.measures = TRUE, standardize = TRUE, rsq = TRUE)
fitMeasures(fit_crk, "df")

# correlation SEM

dat_sem_cor <- cbind(pca_sem_df, microb_axes)
colnames(dat_sem_cor)[19:20] <- c("Micro_Axis.1", "Micro_Axis.2")
cor_sem_model <- '
Chem_Axis.1 ~ Env_pH + Env_Temp + Env_2wPrec + Env_Area + LC_DVLO + LC_CRPHAY + LC_FORSHB
Chem_Axis.2 ~ Env_Temp + Env_BgRun + Env_LAI + Env_Area + LC_DVLO + LC_DVHI + LC_FORSHB + Micro_Axis.1
Chem_Axis.3 ~ Env_AnuPrec + Env_LAI + Env_Slope + LC_DVLO + LC_CRPHAY + LC_FORSHB + LC_HERB + LC_WETLND + Micro_Axis.1 + Micro_Axis.2
Micro_Axis.1 ~ Env_pH + Env_Temp + Env_AnuPrec + Env_2wPrec + Env_LAI + Env_Slope + LC_BARE
Micro_Axis.2 ~ Env_AnuPrec + Env_LAI + LC_DVLO + LC_CRPHAY + LC_FORSHB + Chem_Axis.1 + Chem_Axis.3
'

fit_crk_cor <- sem(
  model          = cor_sem_model, 
  sample.cov     = cor(dat_sem_cor),
  sample.nobs    = 46
)

summary(fit_crk_cor, fit.measures = TRUE, standardize = TRUE, rsq = TRUE)
semPaths(fit_crk_cor, whatLabels = "std", style = "ram", layout = "tree")

nldas_raw <- read.csv('creeks_nldas_raw.csv',row.names = 1)

chem_axes_55['ID'] <- substr(row.names(chem_axes_55), 2 ,4)

sd(aggregate(Axis.1~ID, chem_axes_55, max)[,2] - aggregate(Axis.1~ID, chem_axes_55, min)[,2])
sd(aggregate(Axis.2~ID, chem_axes_55, max)[,2] - aggregate(Axis.2~ID, chem_axes_55, min)[,2])
#by(chem_axes_55['Axis.1'],chem_axes_55$ID,min)
